# ScaleR
ScaleR: An R package for robust scaling


Examples: 

Result <- ScaleR(x,y,method='PLS', inter=0.1 plot=TRUE)

Result_Ridge <- ScaleR(x,y,method='Ridge', inter0.1, plot= TRUE, seed=111, k=3) 

Default method is PLS
Default inter is 0.1
Default plot=TRUE
Default k=0
Default seed=1234
