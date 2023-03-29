# pPAP: plant PKS Analysis and Prediction
# Copyright (c) 2016 Yugo Shimizu
# Bioinformatics Center, Institute for Chemical Research, Kyoto University
# pPAP is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or any later version.
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

def calculate_score_hmm_ACX(areafastafile)
  this_file_path = File.expand_path(File.dirname($0))
  hmmfile = "#{this_file_path}/data/R-ACX_a1a3a4.hmm"
  ave_a, sd_a = 1.34672268907563, 0.421248577788624
  ave_c, sd_c = 1.52890756302521, 1.1028807566774
  ave_x, sd_x = 1.03294117647059, 0.518906280553234
  
  defarr = []
  Bio::FlatFile.auto(areafastafile).each_entry do |ent|
    defarr << ent.definition
  end
  
  # HMM scan
  hmmresults = `hmmscan --max --incE 1 #{hmmfile} #{areafastafile}`
  
  # parse HMM result
  i, j = 0, 0
  scoresarr = []
  score_a, score_c, score_x = 0, 0, 0
  hmmresults.split("//\n").each do |hmmresult|
    hmmresult.split("\n").each do |line|
      if /^Scores for complete sequence/ =~ line || (j > 0 && j < 4)
        j += 1
      elsif j >= 4
        if /\d/ =~ line
          linearr = line.split(" ")
          model, score = linearr.at(8), linearr.at(1).to_f / 17 # divide by the sequence length (=17)
          score = 0 if score < 0
          if model == "R-4-A_a1a3a4"
            score_a = (score - ave_a) / sd_a
          elsif model == "R-4-C_a1a3a4"
            score_c = (score - ave_c) / sd_c
          elsif model == "R-2-X_a1a3a4"
            score_x = (score - ave_x) / sd_x
          end
        else
          scoresarr << [defarr.at(i), score_a, score_c, score_x]
          j = 0
          score_a, score_c, score_x = 0, 0, 0
          break
        end
        j += 1
      end
    end
    i += 1
  end
  return scoresarr
end
