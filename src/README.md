# Source code for HK-IDLISOAnalysis

Notes for development:

- Logic:
  1) Class DataReader should be responsible for reading the input data files into Run/Source objects and PMT objects
  2) Class DataRunInfo should contain all the run-related/source-related info
  3) (Future) Class DataPMT should contain PMT objects with the corresponding PMT-related info
  4) Class DataFilter should apply PMT/data selection/cuts
  5) Class OpticalModel deals with the optical models for analyses, chi2/likelihood functions, etc
  6) (Future) Class Fit should deal with the minimisation of the model parameters (can have different minimisation methods)
  7) (Future) Some class to "store" the results
