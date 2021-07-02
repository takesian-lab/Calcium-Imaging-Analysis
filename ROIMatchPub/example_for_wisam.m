%Example data for wisam

%Let's say I have 4 fall.mats called A, B, C, and D
falls = {'A', 'B', 'C', 'D'};

%Lists of all cells in each fall.mat
A = [1 2 3 4 5];
B = [0 2 4 6 8];
C = [1 3 5 7 9];
D = [0 3 4 8 9];

%Each fall got matched with every other fall
matchAB = [2,4; 5,8];
matchAC = [2,7; 4,5];
matchAD = [4,8; 5,9];
matchBC = [2,9; 4,7];
matchBD = [6,4; 8,9];
matchCD = [1,0; 3,3; 5,8];

% N columns = length(falls);
final_match_list = [  1 nan nan nan;...
                      2   4   7 nan;...
                      3 nan nan nan;...
                      4 nan   5   8;...
                      5   8 nan   9;...
                    nan   0 nan nan;...
                    nan   2   9 nan;...
                    nan   6 nan   4;...
                    nan nan   1   0;...
                    nan nan   3   3];
                    
                    
            