from scipy import *
from pylab import *
import math
import scipy
import numpy
from scipy import stats
import os

#Physical parameters and measured variables
kT=4.1
x=8.

###
#Uncomment the value of asymmetry parameter A you'd like to use for the current simulations.
#NuMA measurement, A=-0.48
A= -0.48

#PRC1 measurment, A=0.04
#A= 0.04
###

#NuMA, D=53,700 nm^2/s
k0=840
#PRC1, D=39,900 nm^2/s
#k0=620


#Monte Carlo timing parameters; deltat <= 0.0001
deltat=0.0001
Time = 10
Steps = int(Time/deltat)

#Not used
Position = 0
Force = 0.005
###


#Calculate and return the slope of a Time/Position input
def Slope(Time, Position):
	Slope = scipy.stats.linregress(Time, Position)
	return Slope[0]

#Calculate the new position of the MT binding domain given the instantaneous Force applied to the protein, and the above rates for hopping, asymmetry, etc...
def Hop(Force, Position):
	kplus = math.exp((Force*(0.5*x+A))/kT)*k0
	Pplus = kplus*deltat
	kminus = math.exp(-(Force*(0.5*x-A))/kT)*k0
	Pminus = kminus*deltat
	RANDOM = np.random.random()
	if RANDOM < Pplus:
		Position = Position + 8
	if RANDOM > 1-Pminus:
		Position = Position - 8
	else:
		Position = Position
	return Position

def HopAP(Force, Position):
	kplus = math.exp((Force*(0.5*x-A))/kT)*k0
	Pplus = kplus*deltat
	kminus = math.exp(-(Force*(0.5*x+A))/kT)*k0
	Pminus = kminus*deltat
	RANDOM = np.random.random()
	if RANDOM < Pplus:
		Position = Position + 8
	if RANDOM > 1-Pminus:
		Position = Position - 8
	else:
		Position = Position
	return Position
	
#Simulate the position vs. time behavior for a given number of steps (calculated above from deltat and total simulation time); returns the time series of a given particle's behavior
def SimulateRunSingleHead(Force, Steps):
	POS=np.empty(Steps)
	POS.fill(0)
	TIME=np.empty(Steps)
	TIME.fill(0)
	for i in range(Steps):
		POS[i]=Hop(Force, POS[i-1])
		TIME[i]=deltat*i
	return TIME,POS



#Sets dimer 'stiffness'; this can vary from ~0.01 (loose) to 0.1 or even 1. (very stiff)
kNuMA = 0.01

###Actual measured values:
#kPRC1 = 0.04
#kNuMA = -0.48
#kEB1 = 0.39

def SimulateRunDimer(Steps,DrivingFunction,kNuMA):
	POS=np.empty((Steps,2))
	POS.fill(0)
	TIME=np.empty(Steps)
	TIME.fill(0)
	ForceT=np.empty(Steps)
	ForceT.fill(0.)
	for i in range(Steps):
		ForceTop=kNuMA*(POS[i-1,1]-POS[i-1,0])
		ForceBottom=kNuMA*(POS[i-1,0]-POS[i-1,1])
		POS[i,0]=Hop(ForceTop, POS[i-1,0])+DrivingFunction[i-1]
		POS[i,1]=Hop(ForceBottom, POS[i-1,1])
		TIME[i]=deltat*i
		ForceT[i]=ForceTop
	return TIME,POS,ForceT


#FLIP MT to antiparallel configuration.	
def SimulateRunDimerAP(Steps,DrivingFunction,kNuMA):
	POS=np.empty((Steps,2))
	POS.fill(0)
	TIME=np.empty(Steps)
	TIME.fill(0)
	ForceT=np.empty(Steps)
	ForceT.fill(0.)
	for i in range(Steps):
		ForceTop=kNuMA*(POS[i-1,1]-POS[i-1,0])
		ForceBottom=kNuMA*(POS[i-1,0]-POS[i-1,1])
		POS[i,0]=Hop(ForceTop, POS[i-1,0])+DrivingFunction[i-1]
		POS[i,1]=HopAP(ForceBottom, POS[i-1,1])
		TIME[i]=deltat*i
		ForceT[i]=ForceTop
	return TIME,POS,ForceT
	
def MultipleDimerTrials(Iterations):
	VEL=np.empty(Iterations)
	VEL.fill(0.)
	for i in range(Iterations):
		TIME,POS = SimulateRunDimer(Steps)
		VEL[i] = Slope(TIME,POS[:,0])
	return VEL

#Simulate the Force vs Velocity curve for a single MT binding domain. Runs from -ForceRange to +ForceRange in increments of ForceRange/ForceStepsTotal (i.e., if ForceRange is from -1 to 1, and ForceStepsTotal=21, it'll choose forces of (-1.0, -0.9, -0.8, ... 0.8, 0.9, 1.0). At each force step, run SimulateRunSingleHead and output the velocity (slope of position time series); repeat this Iter times and take the mean of the velocities; return this and the force value. Whole function returns array of Force and calculated Velocity values.
def SINGLEMTHEAD(ForceStepsTotal,Steps,Iter,ForceRange):
	F=np.empty(ForceStepsTotal)
	F.fill(0.)
	V=np.empty(ForceStepsTotal)
	V.fill(0.)
	for i in range(ForceStepsTotal):
		SLPS=np.empty(Iter)
		SLPS.fill(0)
		Force = ForceRange*(2*(float(i)/float(ForceStepsTotal-1))-1.)
		for j in range(Iter):
			TIME,POS = SimulateRunSingleHead(Force,Steps)
			SLPS[j]=Slope(TIME,POS)
		F[i]=Force
		V[i]=mean(SLPS)
		figure(0)
		plot(TIME,POS,'-')
	return V,F

#Example: Return Vel and Force for Force points (-0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1); calculate two (2) slopes per force point, find mean.
#V,F = SINGLEMTHEAD(11,Steps,2,0.1)

figure(0)
plt.xlabel('Time (seconds)')
plt.ylabel('Distance moved (nm)')

#Uncomment this to run Single MThead simulation for multiple Force values, multiple iterations per force
#V,F = SINGLEMTHEAD(11,Steps,9,0.2)

SineAmp = 8
SineFreq = 50

def GenerateSineWave(Steps,SineAmp,SineFreq):
	SINEWAVE = np.empty(Steps)
	SINEWAVE.fill(0.)
	for i in range(Steps):
		SINEWAVE[i]=(SineAmp*(2*pi*deltat*SineFreq))*sin(2*pi*SineFreq*deltat*i)
	return SINEWAVE

	
def GenerateRandomNoise(Steps,Amp):
	RANDOMNOISE = np.empty(Steps)
	RANDOMNOISE.fill(0.)
	for i in range(Steps):
		RANDOMNOISE[i]=standard_normal()
	RANDOMNOISE=Amp*RANDOMNOISE
	return RANDOMNOISE


ZEROES=np.empty(Steps)
ZEROES.fill(0.)

#Turn this on to run one single simulation of a dimer's response to sinusoidal driving. Sine wave amplitude and frequency can be set in line 1, and the simulation runs with that function, plus a variable "dimer stiffness" term. Two figures are generated; the positions of each of the NuMA MT binding domains, and the force experienced by the top unit of the dimer.
	
#SINEWAVE=GenerateSineWave(Steps,SineAmp,SineFreq)
#TIMEDIMER,POSDIMER,FORCETOP = SimulateRunDimer(Steps,SINEWAVE,0.005)
#figure(0)
#plot(TIMEDIMER,POSDIMER[:,0],'r-')
#plot(TIMEDIMER,POSDIMER[:,1],'b-')
#figure(1)
#plot(TIMEDIMER,FORCETOP,'go')


def MultipleRunsDimerWithNoise(Iterations,Steps,kNuMA,Amp,PlotQuery,VelPlotQuery):
	VEL=np.empty(Iterations)
	VEL.fill(0.)
	for i in range(Iterations):
		RANDOMNOISE=GenerateRandomNoise(Steps,Amp)
		TIMEDIMER,POSDIMER,FORCETOP = SimulateRunDimer(Steps,RANDOMNOISE,kNuMA)
		VEL[i]=Slope(TIMEDIMER,((POSDIMER[:,0]+POSDIMER[:,1])/2))
		DIMERCENTER=(POSDIMER[:,0]+POSDIMER[:,1])/2
		if PlotQuery == 1:
			figure(0)
			plot(TIMEDIMER,DIMERCENTER)		
	if VelPlotQuery == 1:
		figure(1)
		plot(VEL,'o-')
	return VEL


def MultipleRunsDimerWithOscillations(Iterations,Steps,kNuMA,SineAmp,SineFreq,PlotQuery,VelPlotQuery):
	SINEFUNCTION=GenerateSineWave(Steps,SineAmp,SineFreq)
	VEL=np.empty(Iterations)
	VEL.fill(0.)
	for i in range(Iterations):
		TIMEDIMER,POSDIMER,FORCETOP = SimulateRunDimer(Steps,SINEFUNCTION,kNuMA)
		VEL[i]=Slope(TIMEDIMER,((POSDIMER[:,0]+POSDIMER[:,1])/2))
		DIMERCENTER=(POSDIMER[:,0]+POSDIMER[:,1])/2
		if PlotQuery == 1:
			figure(0)
			plot(TIMEDIMER,DIMERCENTER)		
	if VelPlotQuery == 1:
		figure(1)
		plot(VEL,'o-')
	return VEL

def MultipleRunsDimerWithOscillationsAP(Iterations,Steps,kNuMA,SineAmp,SineFreq,PlotQuery,VelPlotQuery):
	SINEFUNCTION=GenerateSineWave(Steps,SineAmp,SineFreq)
	VEL=np.empty(Iterations)
	VEL.fill(0.)
	for i in range(Iterations):
		TIMEDIMER,POSDIMER,FORCETOP = SimulateRunDimerAP(Steps,SINEFUNCTION,kNuMA)
		VEL[i]=Slope(TIMEDIMER,((POSDIMER[:,0]+POSDIMER[:,1])/2))
		DIMERCENTER=(POSDIMER[:,0]+POSDIMER[:,1])/2
		if PlotQuery == 1:
			figure(0)
			plot(TIMEDIMER,DIMERCENTER)		
	if VelPlotQuery == 1:
		figure(1)
		plot(VEL,'o-')
	return VEL

#Determine dimer velocity vs. sine wave amplitude relationship, parallel configuration:	
def FindMeanVelocityVersusSineWaveAmplitude(TotalAmpSteps,Iterations,Steps,kNuMA,SineAmpIncrement,SineFreq):
	directory = os.getcwd()
	filename = directory+'\\SineWave oscillation kNuMA-'+str(kNuMA)+' Freq-'+str(SineFreq)+' N_per-'+str(Iterations)+' Time-'+str(Time)+' deltat-'+str(deltat)+'.dat'
	MEANVEL=np.empty(TotalAmpSteps)
	MEANVEL.fill(0.)
	STDVEL=np.empty(TotalAmpSteps)
	STDVEL.fill(0.)
	SINEWAVEAMPARRAY=np.empty(TotalAmpSteps)
	SINEWAVEAMPARRAY.fill(0.)
	for i in range(TotalAmpSteps):
		VEL = MultipleRunsDimerWithOscillations(Iterations,Steps,kNuMA,(SineAmpIncrement*i),SineFreq,0,0)
		MEANVEL[i] = mean(VEL)
		STDVEL[i] = std(VEL)
		SINEWAVEAMPARRAY[i] = (SineAmpIncrement*float(i))
	PARAMETERS=[kNuMA,SineFreq,Iterations,Time,deltat,Steps]
	PARAMETERS=numpy.asarray(PARAMETERS)
	ALLDATA=[SINEWAVEAMPARRAY,MEANVEL,STDVEL]
	ALLDATA=np.array(ALLDATA)
	ALLDATA=transpose(ALLDATA)
	savetxt(filename,ALLDATA,fmt='%1.5f',delimiter='\t')
	return ALLDATA

	
#Determine dimer velocity vs. sine wave amplitude relationship, ANTIparallel configuration:
def FindMeanVelocityVersusSineWaveAmplitudeAP(TotalAmpSteps,Iterations,Steps,kNuMA,SineAmpIncrement,SineFreq):
	directory = os.getcwd()
	filename = directory+'\\SineWave oscillation AP kNuMA-'+str(kNuMA)+' Freq-'+str(SineFreq)+' N_per-'+str(Iterations)+' Time-'+str(Time)+' deltat-'+str(deltat)+'.dat'
	MEANVEL=np.empty(TotalAmpSteps)
	MEANVEL.fill(0.)
	STDVEL=np.empty(TotalAmpSteps)
	STDVEL.fill(0.)
	SINEWAVEAMPARRAY=np.empty(TotalAmpSteps)
	SINEWAVEAMPARRAY.fill(0.)
	for i in range(TotalAmpSteps):
		VEL = MultipleRunsDimerWithOscillationsAP(Iterations,Steps,kNuMA,(SineAmpIncrement*i),SineFreq,0,0)
		MEANVEL[i] = mean(VEL)
		STDVEL[i] = std(VEL)
		SINEWAVEAMPARRAY[i] = (SineAmpIncrement*float(i))
	PARAMETERS=[kNuMA,SineFreq,Iterations,Time,deltat,Steps]
	PARAMETERS=numpy.asarray(PARAMETERS)
	ALLDATA=[SINEWAVEAMPARRAY,MEANVEL,STDVEL]
	ALLDATA=np.array(ALLDATA)
	ALLDATA=transpose(ALLDATA)
	savetxt(filename,ALLDATA,fmt='%1.5f',delimiter='\t')
	return ALLDATA
	
	
def FindMeanVelocityVersusNoiseAmplitude(TotalAmpSteps,Iterations,Steps,kNuMA,AmpIncrement):
	directory = os.getcwd()
	filename = directory+'\\Noise kNuMA-'+str(kNuMA)+' N_per-'+str(Iterations)+' Time-'+str(Time)+' deltat-'+str(deltat)+'.dat'
	MEANVEL=np.empty(TotalAmpSteps)
	MEANVEL.fill(0.)
	STDVEL=np.empty(TotalAmpSteps)
	STDVEL.fill(0.)
	NOISEAMPARRAY=np.empty(TotalAmpSteps)
	NOISEAMPARRAY.fill(0.)
	for i in range(TotalAmpSteps):
		VEL = MultipleRunsDimerWithNoise(Iterations,Steps,kNuMA,(AmpIncrement*i),0,0)
		MEANVEL[i] = mean(VEL)
		STDVEL[i] = std(VEL)
		NOISEAMPARRAY[i] = (AmpIncrement*float(i))
	PARAMETERS=[kNuMA,Iterations,Time,deltat,Steps]
	PARAMETERS=numpy.asarray(PARAMETERS)
	ALLDATA=[NOISEAMPARRAY,MEANVEL,STDVEL]
	ALLDATA=np.array(ALLDATA)
	ALLDATA=transpose(ALLDATA)
	savetxt(filename,ALLDATA,fmt='%1.5f',delimiter='\t')
	return ALLDATA
	
	
	
def ProduceSingleTraceForPlotting(Steps,SineAmp,SineFreq,kNuMA):
	SINEWAVE=GenerateSineWave(Steps,SineAmp,SineFreq)
	TIME,POS,ForceT=SimulateRunDimer(Steps,SINEWAVE,kNuMA)
	ALLDATA=[TIME,POS[:,0],POS[:,1],ForceT]
	ALLDATA=transpose(ALLDATA)
	directory = os.getcwd()
	filename = directory+'\\SampleTrace kNuMA-'+str(kNuMA)+' Amp_nm-'+str(SineAmp)+' Freq-'+str(SineFreq)+' Time-'+str(Time)+' deltat-'+str(deltat)+'.dat'
	savetxt(filename,ALLDATA,fmt='%1.5f',delimiter='\t')
	plot(TIME,POS)
	return ALLDATA


	
#directory = os.getcwd()
#filename = directory+'\\SineWave oscillation kNuMA-'+str(kNuMA)+' Freq-'+str(SineFreq)+' N_per-'+str(Iterations)+'.dat'



# figure(1)
# plot(V,F,'ro-')
# plt.xlabel('Velocity (nm/s)')
# plt.ylabel('Force (pN)')


show()