% scripts to test block.m
fclose all;
clear

testblock = block(100);
gendata; % generate test data

% as a test, just use default values
testblock.updateMpara;
testblock.updatePpara;
testblock.updateSpara;
% testblock.updateLpara;  

testblock.computeConstants;

testblock.maxRecordRows = 5;
Te = testblock.ppara.Te;
testblock.addElements(10,x,y,z,Te.*ones(1,10));
testblock.updateNBlist;


% needs call the function to end simulation 
% testblock.endsimulation;

