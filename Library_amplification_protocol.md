## Amplification Protocol for High-Throughput Sub-Cloning of Oligo Array Sequences ##  
### Disclaimer this is just what worked for me, there are a lot of ways to skin a cat though :crying_cat_face: ###  
This protocol corresponds to the [Vex-seq Paper](https://doi.org/10.1186/s13059-018-1437-x), and more specifically has details about the library construction shown in A  
![alt text](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5984807/bin/13059_2018_1437_Fig1_HTML.jpg)  
You will need large LB/amp plates (step 8 iii), plus normal stuff for PCR, RT, PCR clean-up, ethanol precipitation, agarose gels, maxi-prep, electrocompetent bacteria, electro-cuvettes, and an electroporator.  
[This paper](https://dx.doi.org/10.1038%2Fnprot.2017.016) has some nice guidelines for a library amplification of GeCKO libraries; similar idea with technical requirements.
1. Order synthetic oligo pool with common primer sequences  
	1. For our project we used Agilent as a supplier, but I've also used CustomArray and I've heard of people using Twist Bioscience as well
1. PCR amplify oligo pool
	It seems important to use lower cycles and I prefer to use GC buffer with Phusion  
	I use a larger reaction to get a higher yield without doing more cycles
	1. 20 uL 5x GC Buffer
	1. 2 uL 10mM dNTPs
	1. 5 uL 10uM FWD Primer
	1. 5 uL 10uM REv Primer
	1. 1 uL Phusion
	1. 2 uL Diluted Oligo Pool (1ng/uL) (I've done this with 5x this concentration and been fine)
	1. 65 uL water  
	98C | 98C 48C 72C | 72C   
	30s | 10s 30s 30s | 5min  
	1x___|__________20x| 1x  
1. Clean up initial PCR reaction with PCR purification columns
1. Digest clean PCR reaction and vector backbone with PstI and XbaI (This is for initial insertion into the vector)
	1. Dephosphorylate vector backbone (usually just add 1 uL of CIP and incubate at 37C for 30 mins)
	1. Run an aliquot of PCR product on 2% agarose gel along with uncut PCR product to ensure cutting was successful
		1. Simultaneously run cut vector on the gel and gel extract the relevant part
1. Ethanol precipitate digested vector and PCR product
	1. Salts can inhibit the ligation reaction
1. Ligate PCR amplified oligo pool into digested vector
	Include an insert(-) control for ligation; I usually do this at a 10-fold dilution from the rest of the reaction
	1. 10 uL 10x NEB ligase Buffer 
	1. 1ug vector (10 uL in this case)
	1. 3:1(insert:vector) molar ratio of insert (5.72 uL of 19.1 ng/uL in this case)
	1. 10 uL T4 ligase
	1. 154.28 uL water
	16C overnight and heat inactivate at 65 for 20 minutes  
1. Ethanol precipitate ligations (resuspend in 21 uL)
	1. This improves downstream transformation efficiency in my hands
1. Electroporate electrocompetent Top10 (7 x 3 uL electroporations)
	It is important to include an insert(-) control here to know background
	1. Right before plating (after recovery) pool all electroporations together
	1. Make 1,000x and 10,000x dilutions of pooled electroporation for separate plating
	1. Plate full concentration on 7 separate large LB/amp plates (I use 150mm petri dishes)
	1. Grow overnight at 37C
1. Count colonies on insert(-) and diluted plates to estimate colony number
	Feng Zhang's paper recommends having 100 colonies/construct, but I've gotten by with 75 before with good library representation
1. Scrape colonies (this gets stinky, so I prefer doing it in a fume hood)
	1. Sterilize a spreader with a flat end to start
	1. Fill up a 50 mL falcon tube with LB and squirt about 5 mLs on a plate
	1. Take the flat end of the spreader and try to scratch off all the colonies so that they are in suspension in the LB
	1. Once most of the colonies are in suspension, tilt the plate and put the LB with the suspended colonies back into the falcon tube
	1. Repeat this for each plate
1. Spin down the falcon tube and do a maxi prep
	I usually get by with one maxi prep/7 plates; if you do more you can run into trouble with having too much debris in the prep and getting a low recovery  
I then repeat this process in principle to insert the exon 3 and part of that intron and insert it in between the test sequence and the barcode (email me if you need more details)  



	
		
	
