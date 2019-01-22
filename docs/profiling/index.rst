Community profiling
===================
The community composition of your samples can be analysed using `MetaPhlAn2 <https://bitbucket.org/biobakery/metaphlan2>`_
(`Truong et al. 2015 <https://www.nature.com/articles/nmeth.3589>`_) which uses clade-specific markers genes to estimate
abundances.

MetaPhlAn2
----------
To run MetaPhlAn2 set :code:`metaphlan2:True` in your configfile.

.. Note:: You will have to run snakemake with the :code:`--use-conda` flag in order for metaphlan2 to work.

Output
^^^^^^
MetaPhlAn2 will produce krona plots for individual samples under results/metaphlan2/*.html and a merged
krona plot with all samples under results/report/metaphlan2/metaphlan2.krona.html. In addition, a heatmap of abundances
with hierarchical clustering of samples and species is created under results/report/metaphlan2/metaphlan2.heatmap.png.
