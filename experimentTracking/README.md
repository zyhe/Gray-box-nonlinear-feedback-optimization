# Time-Varying Optimization with Input Constraints

## Main contents

`├── data/`: store convergence results  
`├── figure/`: store figures  
`├── params/`  
`│   ├── params_0.mat`: parameter file  
`├── tool_functions/`  
│   `├── dynamics.m`: implement system dynamics  
│   `├── extractTVprob.m`: extract problem data  
│   `├── firstTracking.m`: first-order feedback opt.  
│   `├── generateTVprob.m`: generate time-varying problems  
│   `├── grayBoxTracking.m`: gray-box feedback opt.  
│   `├── initialReg.m`: calculate the initial regret  
│   `├── projection.m`: projection function  
│   `├── stochESTracking.m`: stochastic extremum seeking  
│   `├── zerothTracking.m`: model-free feedback opt.  
`├── closedLoopTracking.m`: analyze the closed-loop performance  
`├── gridSearchParams.m`: explore different parameter combinations  
`├── plotResult.m`: plot convergence results  

## How to Run

- Run `closedLoopTracking`. The results will be stored in the `data` and `figure` folders.  
  If `closedLoopTracking(i)` is called, then `params_{i}.mat` is loaded, where `i = 0, 1, ...`. In this case, run `gridSearchParams.m` to generate various files of parameter combinations.
- `plotResult.m` allows adjusting the layout of figures based on the stored data in the `data` folder.
