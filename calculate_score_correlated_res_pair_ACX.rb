# pPAP: plant PKS Analysis and Prediction
# Copyright (c) 2016 Yugo Shimizu
# Bioinformatics Center, Institute for Chemical Research, Kyoto University
# pPAP is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or any later version.
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

#### BEGIN OF SUBFUNCTIONS ####

# align test-set sequence with existing alignment and create position conversion table
def make_test_pos_conv_table(def_seqarr, testsetseq, alnfile, newfile)
  fasta_aln = ""
  def_seqarr.each do |def_seq|
    definition, seq = def_seq
    fasta_aln << seq.to_fasta(definition, 60)
  end
  open(alnfile, "w"){|wf| wf.puts fasta_aln}
  open(newfile, "w"){|wf| wf.puts testsetseq.to_fasta("test_set_seq", 60)}
  def_seqarr_aln = []
  `mafft --quiet --mapout --add #{newfile} --keeplength #{alnfile}`
  pos_conv_table = Hash.new()
  open(newfile + ".map") do |f|
    f.gets
    f.gets
    while line = f.gets
      res, pos_test, pos_aln = line.split(",").collect{|i| i.strip}
      pos_conv_table[pos_aln.to_i - 1] = pos_test.to_i - 1 if /\d/ =~ pos_aln # align -> unalign
    end
  end
  return pos_conv_table
end

# count the number of residue pair in each position pair
def get_residue_pair_count(seqarr, corpos_ijarr, alpha, beta)
  beta = (alpha ** 2) / (seqarr.length + 2 * alpha * 21) unless beta # smoothing parameter for pair
  residue_pair_count = Hash.new()
  corpos_ijarr.each do |ij|
    i, j = ij
    pos_i_arr = seqarr.collect{|seq| seq[i, 1]}
    pos_j_arr = seqarr.collect{|seq| seq[j, 1]}
    pos_ijarr = [pos_i_arr, pos_j_arr].transpose
    residue_pair_count_ij = Hash.new(beta) # add pseuodcount = beta
    pos_ijarr.each{|residue_pair_ij| residue_pair_count_ij[residue_pair_ij] += 1}
    residue_pair_count[ij] = residue_pair_count_ij
  end
  return residue_pair_count
end

# scoring for correlation
def calculate_score_cor(corpos_ijarr, type_residue_pair_count, type_residue_count, type_seqN, testsetseq, alpha, beta, pos_conv_table)
  beta = (alpha ** 2) / (type_seqN + 2 * alpha * 21) unless beta # smoothing parameter for pair
  score = 0
  corpos_ijarr.each do |ij|
    i, j = ij
    if pos_conv_table
      if pos_conv_table[i] && pos_conv_table[j]
        testset_residue_pos_ij = [testsetseq[pos_conv_table[i], 1], testsetseq[pos_conv_table[j], 1]]
      elsif pos_conv_table[i]
        testset_residue_pos_ij = [testsetseq[pos_conv_table[i], 1], "-"]
      elsif pos_conv_table[j]
        testset_residue_pos_ij = ["-", testsetseq[pos_conv_table[j], 1]]
      else
        testset_residue_pos_ij = ["-", "-"]
      end
    else
      testset_residue_pos_ij = [testsetseq[i, 1], testsetseq[j, 1]]
    end
    testset_residue_pos_i, testset_residue_pos_j = testset_residue_pos_ij
    pair_freq_ij = type_residue_pair_count[ij][testset_residue_pos_ij] / (type_seqN + 441 * beta) # Additive smoothing, 441 = 21*21 (#21 = 20aa + gap)
    freq_i = type_residue_count[i][testset_residue_pos_i] / (type_seqN + 21 * alpha) # Additive smoothing (#21 = 20aa + gap)
    freq_j = type_residue_count[j][testset_residue_pos_j] / (type_seqN + 21 * alpha) # Additive smoothing (#21 = 20aa + gap)
    score += Math.log2(pair_freq_ij / (freq_i * freq_j))
  end
  score /= corpos_ijarr.length
  return score
end

#### END OF SUBFUNCTIONS ####

#### MAIN FUNCTION ####
def calculate_score_cor_ACX(queryfastaarr, tempdir)
  #### BEGIN OF PARAMETER SETTINGS ####
  this_file_path = File.expand_path(File.dirname($0))
  alnfile, newfile = "#{tempdir}/pks3_aln.fasta", "#{tempdir}/pks3_new.fasta" # avoid collision of filename in parallel execution
  modelfiles = ["#{this_file_path}/data/R-4-A_aln.fasta", "#{this_file_path}/data/R-4-C_aln.fasta", "#{this_file_path}/data/R-2-X_aln.fasta"] # reference FASTA alignment file for each reaction type
  corposfiles = ["#{this_file_path}/data/R-4-A_corpair.txt", "#{this_file_path}/data/R-4-C_corpair.txt", "#{this_file_path}/data/R-2-X_corpair.txt"] # correlated positions file
  ave = [-2.61144729899128, -2.76589957977704, -3.41939061831597] # for normalization
  sd = [1.52612586504038, 1.88324180308975, 1.39674257969728] # for normalization
  alpha = 1.0 / 21.0 # Smoothing parameter
  beta = nil # Smoothing parameter, nil -> calculated from alpha
  #### END OF PARAMETER SETTINGS ####

  #### BEGIN OF PROCESS ####
  # load query data
  def_seqarrq = []
  queryfastaarr.each do |queryfasta|
    ent = Bio::FastaFormat.new(queryfasta)
    def_seqarrq << [ent.definition, ent.seq]
  end
  
  model_scorearr = Hash.new{|h, k| h[k] = []}
  modelfiles.each_with_index do |modelfile, m|
    def_seqarr_aln_c = []
    Bio::FlatFile.auto(modelfile).each_entry{|ent| def_seqarr_aln_c << [ent.definition, ent.seq]}
    
    # read correlated position-pairs
    corpos_ijarr_aln_c = []
    open(corposfiles.at(m)) do |f|
      while line = f.gets
        corpos_ijarr_aln_c << line.chomp.split("\t").collect{|x| x.to_i}
      end
    end
    
    seqarr_c = def_seqarr_aln_c.collect{|def_seq| def_seq.at(1)}
    residue_pair_count_c = get_residue_pair_count(seqarr_c, corpos_ijarr_aln_c, alpha, beta)
    seqnum_c = seqarr_c.length
    
    residue_count_c = Hash.new()
    corpos_ijarr_aln_c.flatten.uniq.each do |pos|
      residue_count_c_i = Hash.new(alpha) # add pesudocount = alpha
      seqarr_c.collect{|seq| seq[pos, 1]}.each{|res| residue_count_c_i[res] += 1}
      residue_count_c[pos] = residue_count_c_i
    end
    
    def_seqarrq.each do |def_seq|
      querydef, queryseq = def_seq
      # calculate score of query for model type
      pos_conv_table = make_test_pos_conv_table(def_seqarr_aln_c, queryseq, alnfile, newfile)
      score_c = calculate_score_cor(corpos_ijarr_aln_c, residue_pair_count_c, residue_count_c, seqnum_c, queryseq, alpha, beta, pos_conv_table)
      model_scorearr[m] << (score_c - ave[m]) / sd[m]
    end
  end
  
  scoreshash = Hash.new()
  def_seqarrq.each_with_index do |def_seq, i|
    scoreshash[def_seq.at(0)] = [model_scorearr[0].at(i), model_scorearr[1].at(i), model_scorearr[2].at(i)]
  end
  return scoreshash
end
