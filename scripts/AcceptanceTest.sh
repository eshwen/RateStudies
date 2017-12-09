### tests VBF seed on VBFHH2B2Taus signal wrt DoubleIsoTau31e
./AcceptanceTest.py 25 70 --ditaupt 32 --VBFHH --VBFsubL1 35 --VBFleadL1 110 --VBFMjjL1 620

# tests on VBF seeds @ 2E34 wrt DoubleIsoTau31er (emu)
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 30 --VBFleadL1 90 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 30 --VBFleadL1 100 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 35 --VBFleadL1 100 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 40 --VBFleadL1 115 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 40 --VBFleadL1 110 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 45 --VBFleadL1 110 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 45 --VBFleadL1 115 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 50 --VBFleadL1 115 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 40 --VBFleadL1 120 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 45 --VBFleadL1 120 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 35 --VBFleadL1 110 --VBFMjjL1 630
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 35 --VBFleadL1 115 --VBFMjjL1 630
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 35 --VBFleadL1 120 --VBFMjjL1 630
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 35 --VBFleadL1 110 --VBFMjjL1 640
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 35 --VBFleadL1 115 --VBFMjjL1 640
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 35 --VBFleadL1 120 --VBFMjjL1 640
#for rej
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 35 --VBFleadL1 110 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --emu --VBFsubL1 30 --VBFleadL1 115 --VBFMjjL1 620


# tests on VBF seeds @ 2E34 wrt DoubleIsoTau31er (unpack)
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --VBFsubL1 30 --VBFleadL1 90 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --VBFsubL1 30 --VBFleadL1 100 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --VBFsubL1 35 --VBFleadL1 100 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --VBFsubL1 40 --VBFleadL1 115 --VBFMjjL1 620
#for rej
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --VBFsubL1 35 --VBFleadL1 110 --VBFMjjL1 620
#./AcceptanceTest.py 25 70 --ditaupt 31 --VBF --VBFsubL1 30 --VBFleadL1 115 --VBFMjjL1 620           
### seeds 14kHz at 18E34 wrt DoubleIsoTau30er
#VBF signal
# ./AcceptanceTest.py 25 70 --offptpair 80 --VBF --VBFtag
# ./AcceptanceTest.py 25 43 --VBF --VBFtag
# ./AcceptanceTest.py 26 33 --VBF --VBFtag
# ./AcceptanceTest.py 27 23 --VBF --VBFtag
# ./AcceptanceTest.py 28 15 --VBF --VBFtag
# ./AcceptanceTest.py 28 15 --ggH 

#ggH signal		 
# ./AcceptanceTest.py 25 70 --offptpair 80 --ggH --boosted
# ./AcceptanceTest.py 25 43 --ggH --boosted
# ./AcceptanceTest.py 26 33 --ggH --boosted
# ./AcceptanceTest.py 27 23 --ggH --boosted
# ./AcceptanceTest.py 28 15 --ggH --boosted
			 
#Ztt signal		 
# ./AcceptanceTest.py 25 70 --offptpair 80 --Ztt
# ./AcceptanceTest.py 25 43 --Ztt 
# ./AcceptanceTest.py 26 33 --Ztt 
# ./AcceptanceTest.py 27 23 --Ztt 
# ./AcceptanceTest.py 28 15 --Ztt 
			 
#HHbbtt
# ./AcceptanceTest.py 25 70 --offptpair 80 --HHbbtt
# ./AcceptanceTest.py 25 43 --HHbbtt 
# ./AcceptanceTest.py 26 33 --HHbbtt 
# ./AcceptanceTest.py 27 23 --HHbbtt 
# ./AcceptanceTest.py 28 15 --HHbbtt

### seeds 14kHz at 2E34 wrt DoubleIsoTau31er
#VBF signal
# ./AcceptanceTest.py 25 70 --ditaupt 31 --offptpair 80 --VBF --VBFtag
# ./AcceptanceTest.py 25 49 --ditaupt 31 --VBF --VBFtag
# ./AcceptanceTest.py 26 41 --ditaupt 31 --VBF --VBFtag
# ./AcceptanceTest.py 27 31 --ditaupt 31 --VBF --VBFtag
# ./AcceptanceTest.py 28 22 --ditaupt 31 --VBF --VBFtag

#ggH signal
# ./AcceptanceTest.py 25 70 --ditaupt 31 --offptpair 80 --ggH --boosted
# ./AcceptanceTest.py 25 49 --ditaupt 31 --ggH --boosted
# ./AcceptanceTest.py 26 41 --ditaupt 31 --ggH --boosted
# ./AcceptanceTest.py 27 31 --ditaupt 31 --ggH --boosted
# ./AcceptanceTest.py 28 22 --ditaupt 31 --ggH --boosted

#Ztt signal
# ./AcceptanceTest.py 25 70 --ditaupt 31 --offptpair 80 --Ztt
# ./AcceptanceTest.py 25 49 --ditaupt 31 --Ztt 
# ./AcceptanceTest.py 26 41 --ditaupt 31 --Ztt 
# ./AcceptanceTest.py 27 31 --ditaupt 31 --Ztt 
# ./AcceptanceTest.py 28 22 --ditaupt 31 --Ztt 

#HHbbtt
# ./AcceptanceTest.py 25 70 --ditaupt 31 --offptpair 80 --HHbbtt
# ./AcceptanceTest.py 25 49 --ditaupt 31 --HHbbtt 
# ./AcceptanceTest.py 26 41 --ditaupt 31 --HHbbtt 
# ./AcceptanceTest.py 27 31 --ditaupt 31 --HHbbtt 
# ./AcceptanceTest.py 28 22 --ditaupt 31 --HHbbtt

### seeds 14kHz at 2E34 wrt DoubleIsoTau28er
#VBF signal
# ./AcceptanceTest.py 25 70 --ditaupt 28 --offptpair 80 --VBF --VBFtag
# ./AcceptanceTest.py 25 49 --ditaupt 28 --VBF --VBFtag
# ./AcceptanceTest.py 26 41 --ditaupt 28 --VBF --VBFtag
# ./AcceptanceTest.py 27 31 --ditaupt 28 --VBF --VBFtag
# ./AcceptanceTest.py 28 22 --ditaupt 28 --VBF --VBFtag

#ggH signal
# ./AcceptanceTest.py 25 70 --ditaupt 28 --offptpair 80 --ggH --boosted
# ./AcceptanceTest.py 25 49 --ditaupt 28 --ggH --boosted
# ./AcceptanceTest.py 26 41 --ditaupt 28 --ggH --boosted
# ./AcceptanceTest.py 27 31 --ditaupt 28 --ggH --boosted
# ./AcceptanceTest.py 28 22 --ditaupt 28 --ggH --boosted

#Ztt signal
# ./AcceptanceTest.py 25 70 --ditaupt 28 --offptpair 80 --Ztt
# ./AcceptanceTest.py 25 49 --ditaupt 28 --Ztt 
# ./AcceptanceTest.py 26 41 --ditaupt 28 --Ztt 
# ./AcceptanceTest.py 27 31 --ditaupt 28 --Ztt 
# ./AcceptanceTest.py 28 22 --ditaupt 28 --Ztt 

#HHbbtt
# ./AcceptanceTest.py 25 70 --ditaupt 28 --offptpair 80 --HHbbtt
# ./AcceptanceTest.py 25 49 --ditaupt 28 --HHbbtt 
#./AcceptanceTest.py 26 41 --ditaupt 28 --HHbbtt 
# ./AcceptanceTest.py 27 31 --ditaupt 28 --HHbbtt 
# ./AcceptanceTest.py 28 22 --ditaupt 28 --HHbbtt


### HighBoost or DiTau wrt DoubleIsoTau29er

# ./AcceptanceTest.py 25 70 --ditauptOR 31 --ditaupt 29 --offptpair 80 --VBF --VBFtag
# ./AcceptanceTest.py 25 70 --ditauptOR 31 --ditaupt 29 --offptpair 80 --ggH --boosted
# ./AcceptanceTest.py 25 70 --ditauptOR 31 --ditaupt 29 --offptpair 80 --HHbbtt
# ./AcceptanceTest.py 25 70 --ditauptOR 31 --ditaupt 29 --offptpair 80 --Ztt



### High Boost Seeds seeds TOTAL 14kHz at 2E34 wrt DoubleIsoTau31er
#VBF signal
# ./AcceptanceTest.py 25 83 --ditauptOR 34 --ditaupt 31 --VBF --VBFtag
# ./AcceptanceTest.py 26 87 --ditauptOR 33 --ditaupt 31 --VBF --VBFtag
# ./AcceptanceTest.py 27 79 --ditauptOR 33 --ditaupt 31 --VBF --VBFtag
# ./AcceptanceTest.py 28 84 --ditauptOR 32 --ditaupt 31 --VBF --VBFtag

#ggH signal

# ./AcceptanceTest.py 25 83 --ditauptOR 34 --ditaupt 31 --ggH --boosted
# ./AcceptanceTest.py 26 87 --ditauptOR 33 --ditaupt 31 --ggH --boosted
# ./AcceptanceTest.py 27 79 --ditauptOR 33 --ditaupt 31 --ggH --boosted
# ./AcceptanceTest.py 28 84 --ditauptOR 32 --ditaupt 31 --ggH --boosted

#Ztt signal

# ./AcceptanceTest.py 25 83 --ditauptOR 34 --ditaupt 31 --Ztt 
# ./AcceptanceTest.py 26 87 --ditauptOR 33 --ditaupt 31 --Ztt 
# ./AcceptanceTest.py 27 79 --ditauptOR 33 --ditaupt 31 --Ztt 
# ./AcceptanceTest.py 28 84 --ditauptOR 32 --ditaupt 31 --Ztt 
						        
#HHbbtt						        
						        
# ./AcceptanceTest.py 25 83 --ditauptOR 34 --ditaupt 31 --HHbbtt 
# ./AcceptanceTest.py 26 87 --ditauptOR 33 --ditaupt 31 --HHbbtt
# ./AcceptanceTest.py 27 79 --ditauptOR 33 --ditaupt 31 --HHbbtt 
# ./AcceptanceTest.py 28 84 --ditauptOR 32 --ditaupt 31 --HHbbtt



