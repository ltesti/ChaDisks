# do all the analysis
for i in ./*; do
  if [ -d "$i" ]; then 
    cd "$i"; 
    #python analysis_main.py -d --cpu=20 -s0123456789  # do all analysis
    python analysis_main.py -s8  # prepare results_##.txt
    #python analysis_main.py -s7  # do uv-plot for set of models
    cd ..; 
  fi; 
done


