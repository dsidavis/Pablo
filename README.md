Pablo gave us
  DTL.R
  Farms_to_D.csv
  matrix_to_D.csv

We extracted 
  the functions into functions.R
  the setup for the data inputs into setup.R


model1.R - initial version
model2.R - split stocks on substring() rather than grep
model3.R - avoid substring() and "know" the way to split the stocks.
model4.R - testing if the inputs change across calls


The %*% is the dominant function at 60%.
With nsims = 1, 14067 calls to the model() function.

