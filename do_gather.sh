# gather the analsys results

mkdir ./all_plots
mkdir ./all_bestfit
#mkdir ./chain_diag

for i in ./*; do
  if [ -d "$i" ]; then
    cp -r "$i"/analysis/plots/*.pdf ./all_plots/;
    cp -r "$i"/analysis/bestfit/bestfit_*.txt ./all_bestfit/;
    cp "$i"/analysis/result*.txt ./all_bestfit/;
#    cp -r "$i"/analysis/chain_diag/ ./chain_diag/"$i"/
fi;
done


