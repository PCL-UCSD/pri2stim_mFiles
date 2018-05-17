function myseed = load_rng_seed
% used to ensure same seed used for all randomization analyses. due to the
% way randomization may differ between specific analyses, not guaranteed to
% exactly randomize the same way across analyses, but this at least ensures
% some level of consistency/fairness.

myseed = 83224561; % colleague MS provided over GChat on 9/5/2016


return