# formind-pfts

This repo contains predominantly R code duplicating the FORMIND parameterisation of plant functional types. The main point is to provide quick tools to visualise a PFT with given parameters, including the geometry of the PFT and the growth of the PFT under simplified conditions (constant illumination, no water or temperature constraint).

The R code to generate the PFT geometry and growth is in `FORMIND_tree_models.R`. There are two example folders attempting to match the implementation in R to FORMIND outputs. 

#### Tropical Forest

This uses a modified version of the tropical forest parameterisation provided with FORMIND. It has been switched to use the simpler Chanter parameterisation for diameter increments. Additionally, all mortality and seedlings have been turned off and  the `.pin` file has been edited to start with a single sapling of each of the three PFTs in widely separated patches. This gives a single individual growing without shading or competition and makes it easy to study the growth pattern by looking at the contents of the `.res` file for the simulation.

The R file `tropicalForest_noseeds_chanter.R` goes through a process of:

1. Matching the geometry of the PFTs from the `.res` file to the functions in R.
2. Checking the growth mechanism from the FORMIND `.res` file (that $D_{n+1} = D_{n} + DInc{n+1}$ etc.)
3. Simulating the growth using R and comparing that to the `.res` file. The PDF in the folder shows the differences between the variables calculated in FORMIND and through R are very similar and that the productivity curves through time using the two methods are pretty much in agreement.

There are some unfixed issues - a small but systematic error in the PB calculation and the Chanter model for the intermediate PFT doesn't seem to behave in the same way, although the parameters are the same.

#### Lonesome

I set this up to test a single tree and then got some odd behaviour in the `.res` file, which I've been trying to understand. The R file `lonesome.R` goes through the comparison of the R code to the `.res` data. The key points are:

1. The production curves from the `.res` file showed **zero** maintenance respiration for the first 25 years of growth but growth continues. This then spikes up and growth respiration and growth drop (Top left panel in the `Lonesome_Growth.pdf`).   

    However,`R_Growth`, `R_Main` and `BInc` do still sum to give `PB`. That's what I'd expect - PB is the budget to be split between growth, growth respiration and maintenance respiration.

2. The basic issue is obvious: the growth parameterisation is wrong so  that growth exceeds the available production. If I naively simulate it in R (top right panel in the PDF), negative maintenance respiration is needed to make the budget sum up.

3. I looked at the relationship between `D` and `DInc`, in the `.res` file to see how growth was proceeding and there is a clear truncation of the expected polynomial shape until the tree gets $D=0.6$ at about 25 years.

4. I can reproduce this behaviour pretty closely in R by finding a DInc value that gives growth costs (BInc + R_Growth) that equal the available PB, leaving zero maintenance respiration. If I run this constrained simulation then the production curves match what I see in FORMIND and the simualted `DInc` also matches closely to the values reported in the `.res` file.
