%                     timeStamps: (struct) time stamps for velocity and
%                                 acceleration threshold crossings
% 
%                                 struct fields:
% 
%                                 1) timeStamps.thresholdATimeStamps (double) vector - Dim: number of potential bouts x 2
%                                       A 2 column array where the first column contains upward 
%                                       crossings and the second contains downward crossings
%                                       through threshold A
% 
%                                 2) timeStamps.thresholdBTimeStamps
%                                       (double) vector - Dim: number of potential bouts x 2
%                                       A 2 column array where the first column contains 
%                                       upward crossings and the second contains downward crossings
%                                       through threshold B
% 
%                                 3) timeStamps.nonmotorBouts