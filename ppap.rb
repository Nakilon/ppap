#!/usr/bin/env ruby
# pPAP: plant PKS Analysis and Prediction
# Copyright (c) 2016 Yugo Shimizu
# Bioinformatics Center, Institute for Chemical Research, Kyoto University
# pPAP is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or any later version.
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

this_file_path = File.expand_path(File.dirname($0))
$LOAD_PATH.push(this_file_path)
require "bio"
require "tmpdir"
require "calculate_score_hmm_ACX.rb"
require "calculate_score_correlated_res_pair_ACX.rb"

#### BEGIN OF OPTION SETTINGS ####
require 'optparse'
opt = OptionParser.new
scriptname = File.basename(__FILE__)
Version = "1.1"
usage_text = "Usage:\n% ruby #{scriptname} query.fasta"
opt.banner = usage_text
opt.parse(ARGV)
#### END OF OPTION SETTINGS ####

#### BEGIN OF PARAMETER SETTINGS ####
area_reference_aln_file = "#{this_file_path}/data/R-ACX_aln.fasta"
range = "join(118..122,179..184,291..296)"
errtext = ""
ld_a = [-0.3849737, -1.4807169, 0.4098266, 1.8676184, -0.1400028, -0.9704928]
ld_a_cons = -1.666329
ld_c = [-0.25332139, 0.2200935, -0.12345235, -0.19075879, 1.77633298, 0.03152257]
ld_c_cons = -0.3505823
ld_x = [-0.3279558, 0.5914049, -0.3179594, -0.5722581, -0.8447108, 3.3680182]
ld_x_cons = -3.758074
#### END OF PARAMETER SETTINGS ####

#### BIGIN OF PROCESS ####
if ARGV.length == 0
  puts usage_text
  exit(1)
end
begin
  textinput = ARGF.read

  fastaarr = []
  defarr = []
  definition, seq = nil, ""
  textinput.gsub!(/\r\n|\n|\r/, "\n")
  textinput.split("\n").each do |line|
    if /^>(.+)/ =~ line
      if defarr.include?(definition)
        errtext << "# Caution! Redundant definition: #{definition.slice(2..-1)}" + "\n"
      else
        if definition # ignore lines before first definition line
          fastaarr << Bio::Sequence::AA.new(seq).to_fasta(definition, 60) # input FASTA of the previous definition
          defarr << definition
        end
      end
      new_def = "Q_" + $1 # avoid collision of definition
      definition, seq = new_def, ""
    else
      seq << line
    end
  end
  if defarr.include?(definition)
    errtext << "# Caution! Redundant definition: #{definition.slice(2..-1)}" + "\n"
  else
    fastaarr << Bio::Sequence::AA.new(seq).to_fasta(definition, 60)
    defarr << definition
  end

  areafastaarr = []
  prtarr = []
  Dir.mktmpdir(["pks3c", ""]) do |tempdir|
    queryfile = tempdir + "/query.fasta"
    queryareafile = tempdir + "/query_area.fasta"
    open(queryfile, "w"){|f| fastaarr.each{|fasta| f.puts fasta}} # write the query into a fasta file
    Bio::MAFFT::Report.new(`mafft --quiet --add #{queryfile} --keeplength #{area_reference_aln_file}`).entries.each do |ent|
      areafastaarr << ent.seq.splicing(range).to_fasta(ent.definition.slice(2..-1), 60) if defarr.include?(ent.definition)
    end
    open(queryareafile, "w"){|f| areafastaarr.each{|fasta| f.puts fasta}} # write the query into a fasta file
    hmmscoresarr = calculate_score_hmm_ACX(queryareafile)
    corscoreshash = calculate_score_cor_ACX(fastaarr, tempdir)
    hmmscoresarr.each do |definition_hmmscoresACX|
      definition, *hmmscores = definition_hmmscoresACX
      corscores = corscoreshash["Q_" + definition]
      scores = hmmscores + corscores
      lda_score = [ld_a_cons, ld_c_cons, ld_x_cons] # R-4-A, R-4-C, R-2-X
      scores.each_with_index do |score, i|
        lda_score[0] += score * ld_a.at(i)
        lda_score[1] += score * ld_c.at(i)
        lda_score[2] += score * ld_x.at(i)
      end
      if lda_score.at(2) > 0
        prtarr << "R-2-X"
      elsif lda_score.at(0) > 0
        prtarr << "R-4-A"
      elsif lda_score.at(1) > 0
        prtarr << "R-4-C"
      else
        prtarr << "Other"
      end
    end
  end
rescue => e
  puts "Error: Please confirm the input and try again."
  puts e.to_s # debug
  exit(1)
end

# output
print errtext
puts [areafastaarr.collect{|fasta| "Query\t" + fasta[/[^>\n]+\n/]}, prtarr.collect{|prt| "Type\t" + prt}].transpose.collect{|fasta_prt| fasta_prt.join("")}.join("\n//\n")
puts "//"
