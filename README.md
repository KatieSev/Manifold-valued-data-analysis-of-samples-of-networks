# Manifold-valued-data-analysis-of-samples-of-networks
Code to produce results in [1].

The latest code is available at: 

https://github.com/KatieSev/Manifold-valued-data-analysis-of-samples-of-networks

This work was supported by the Engineering and Physical Sciences Research Council [grant number EP/T003928/1 and 
EP/M02315X/1]. The authors are grateful to Michaela Mahlberg, Viola Wiegand and Anthony Hennessey for their help and discussions about the data.
The data are obtained from CLiC https://clic.bham.ac.uk (see [2]). 

[1] Severn, K., Dryden, I.L., and Preston, S.P. (2021). "Manifold valued data analysis of samples of networks, with applications in corpus linguistics", Annals of Applied Statistics (to appear). 

[2] Mahlberg, M., Stockwell, P., De Joode, J., Smith, C. and O’Donnell, M. B. (2016). "CLiC Dickens: novel uses of concordances for the integration of corpus stylistics and cognitive poetics", Corpora 11 433-463.

# Example 

To produce the results in [1], run

source('main_script_to_run.R') 

Only including the top 100 words in the graph Laplcians has been set as the default, for efficiency,
but this may be changed to 1000 (as in the paper) by changing 'topwords<-1000' in main_script_to_run.R.

