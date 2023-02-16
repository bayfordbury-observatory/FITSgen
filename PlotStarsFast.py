import numpy as np
from PIL import Image
import math
from astropy.io import fits

name="m103-b"

#I (e-) = =2.5*POWER(2.5118864315,(21.69885062-J23))

gain=3.947

#read nose
rn = 1.5*4*1.5*4
#dark current
dc = 1

width = 2394
height= 1597


seeing=1.3 #sigma

exp = 15
expin = 1

scale=45*5
offset=25

pedestal = 100-rn

gRange=int(seeing*10)#30

print("gRange " +str(gRange))

starIm = np.zeros((height, width), dtype=float)
	

def main(skymag):





	#res = 100
	
	#sqrtArr = np.zeros((gRange*2+1, gRange*2+1), dtype=np.float)
	
	#for dx in range(-gRange, gRange+1):
	#			for dy in range(-gRange, gRange+1):
	#				
	#				sqrtArr[dx+gRange, dy+gRange] = math.hypot(float(dx), float(dy))



	#skymag=22

	print("skymag " +str(skymag))

	bg=exp*2.5*pow(2.512,(21.16088324-skymag))

	print("bg "+ str(bg))

	im = (starIm*exp/expin )+bg+rn+dc


 
	#print("orig")

	#print(im.shape)
	#print(im)

	print(np.mean(im))
	print(np.std(im))

	#noisy = np.random.poisson(im)

	print("Adding noise")
	
	noisy = np.random.poisson(lam=im, size=None)
	#noisy = im

	#print("noisy")

	#print(noisy.shape)
	#print(noisy)
	print(np.mean(noisy))
	print(np.std(noisy))

	#scaled = np.clip(((noisy-bg)*scale/exp)+offset, 0, 255);

	#print("scaled");


	#print(scaled)

	#eightBit = np.uint8(scaled)
	sixteenbit = np.uint16(np.clip((noisy/gain)+pedestal, 0, 65535))
	
	#print("8 bit");

	#print(eightBit)

	#img = Image.fromarray(eightBit, 'L')

	#img.save('out.png')
		
	hdu = fits.PrimaryHDU(sixteenbit)
	
	#hdu.writeto("stars2_s"+str(skymag)+"_e"+str(exp)+"_sig"+str(seeing)+".fits")
	#hdu.writeto("2mass_sim_s"+str(skymag)+"_e"+str(exp)+"_sig"+str(seeing)+".fits")
	hdu.writeto(name+"_s"+str(skymag)+"_e"+str(exp)+"_sig"+str(seeing)+"_new_subp_fast.fits")
	
#import cProfile, pstats
#profiler = cProfile.Profile()
#profiler.enable()



print("Reading input")

#stars = open("stars.csv", "r")
stars = open(name+".csv", "r")
#stars = open("2mass_v_45.csv", "r")
starLines = stars.readlines()

print("Adding stars")

for line in starLines:
#for line in range(0, 100):

	
	#parts = starLines[line].split(",")
	parts = line.split(",")

	xc = float(parts[0])
	yc = float(parts[1])
	mag = float(parts[2])
	I =2.5*math.pow(2.5118864315,(21.69885062-mag))
	#I = float(parts[2])	
	#mag = float(parts[3])	
	
	#print(I)
	
	
	if xc>-gRange and xc<(width+gRange) and yc>-gRange and yc<(height+gRange):
				
		xs = int(xc-gRange )
		xe = int(xc+gRange)
		ys = int(yc-gRange)
		ye = int(yc+gRange)
		
		if xs<0:
			xs=0
		if ys<0:
			ys=0
		if xe>(width-1):
			xe=width-1
		if ye>(height-1):
			ye=height-1
		
		for x in range(xs, xe):
			for y in range(ys,ye):
				dx = float(x)-xc
				dy = float(y)-yc
				d = math.hypot(dx, dy)#math.sqrt(dx*dx+dy*dy)
				#d = sqrtArr[dx+gRange, dy+gRange]
				#d=1;
				k1 = seeing*2.50662827463
				k2 = -0.5*(d/seeing)*(d/seeing)
				
				#p = (0.143133580102272*I/k1*math.exp(k2))+(0.052276908*0.5*I/k1*math.exp(0.25*k2))
				p = (0.143133580102272*I/k1*np.e**k2)+(0.052276908*0.5*I/k1*np.e**(0.25*k2))
				
				#=1;
				
				#im[x,y]=im[x,y]+p
				starIm[y,x]=starIm[y,x]+p


print("..done")

for skymag in range(15,22):

	main(skymag)
#profiler.disable()
#stats = pstats.Stats(profiler).sort_stats('ncalls')
#stats.print_stats()
