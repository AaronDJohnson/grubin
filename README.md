This repo currently has two functions:
  * Gelman-Rubin Rhat statistic
  * `find_files()`

By using `find_files()`, you can search a subdirectory and find out if the chains output by PTMCMCSampler satisfy the threshold value of 1.1 in the Gelman Rubin statistic.