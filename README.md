# Statistical_Inference_Coursework_Project
This was a teamwork project with 3 of my classmates in 2017. The aim of this project is understanding the process of the Metropolis-Hastings (MH) sampler. 
By using the sampler we investigated the toxicity of a new pharmaceutical treatment which was in the Pahse I clinical trial. 
The sampler allowed us to have a better understanding of the the situationa nd parameters. Only 40 patients were enlisted 
and 16 of them have biomarker. They were randomly split into 10 groups with 4 in each group. Patiets in the same group 
receives the same dose amount (ranging form 0 to 200, with increasement unit around 20). The method for estimating parameters 
is MCMC and the model we used is:
$$
r_i = E_0 + \frac{d^{\lambda}_i E_max}{d^{\lambda}+(ED_50 + \beta x_i)^{\lambda}}+\epsilon_i
$$
