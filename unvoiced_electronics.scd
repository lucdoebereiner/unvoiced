/*

"UNVOICED" FOR VOICE AND LIVE-ELECTRONICS
(c) 2016 by Luc Döbereiner
requires SuperCollider at least 3.6.6. with sc-plugins + PhinGen UGen by Luc Döbereiner

*/

(
o = Server.local.options;
o.numWireBufs = 256;
o.memSize = 262144;
o.blockSize = 16;
s.waitForBoot({

	// adjust path
	~trans = Buffer.read(s, "/home/luc/Work/voice_elec/sec1trans2.wav");
	~input = 0; // adjust to audio interface if necessary
	// first of two adjacent output audio channels
	~output = 2; // adjust to audio interface if necessary
	~outputVc = 0; // output for amplified voice

	// BUFFERS AND BUSSES

	~b1 = Buffer.alloc(s, 2048, 1);
	~b2 = Buffer.alloc(s, 2048, 1);
	~b3 = Buffer.alloc(s, 2048, 1);
	~b4 = Buffer.alloc(s, 2048, 1);
	~b5 = Buffer.alloc(s, 2048, 1);
	~b6 = Buffer.alloc(s, 2048, 1);
	~b7 = Buffer.alloc(s, 2048, 1);
	~b8 = Buffer.alloc(s, 2048, 1);

	~recb1 = Buffer.alloc(s, 44100 * 2.2, 1);
	~recb2 = Buffer.alloc(s, 44100 * 2.2, 1);
	~recb3 = Buffer.alloc(s, 44100 * 2.2, 1);
	~recb4 = Buffer.alloc(s, 44100 * 2.2, 1);

	~wind = Buffer.sendCollection(s, Env([0.0,1,1,0.0],[0.001,1,0.001]).discretize, 1);
	~window = Env.linen(0.02,0.96,0.02,1,[2,1,-2]);
	~winBuf = Buffer.sendCollection(s, ~window.discretize, 1);

	~recWindow = Env.linen(0.02,0.96,0.02,1,[2,1,-2]);
	~recWinBuf = Buffer.sendCollection(s, ~recWindow.discretize, 1);
	~fadeWindow = Env.linen(0.45,0.1,0.45,1,[2,1,-2]);
	~fadeWinBuf = Buffer.sendCollection(s, ~fadeWindow.discretize, 1);

	~voiceBus = Bus.audio(s,1);

	~recS21 = Buffer.alloc(s, 44100 * 2, 1);
	~recS22 = Buffer.alloc(s, 44100 * 2, 1);
	~recS31 = Buffer.alloc(s, 44100 * 4, 1);
	~recS32 = Buffer.alloc(s, 44100 * 4, 1);

	~analBuf1 = Buffer.alloc(s, 2048);
	~analBuf2 = Buffer.alloc(s, 2048);
	~analBuf3 = Buffer.alloc(s, 2048);
	~analBuf4 = Buffer.alloc(s, 2048);

	~recS4 = Buffer.alloc(s, 44100 * 4, 1);

	// SYNTHDEFS

	SynthDef("unvoicedSec1", {arg out=0,inbus=2,bufIn=0,busout,
		rmod=0.05,lagtime=20,amp=0.1,rmin=3,rmax=16,fbAmp=0,
		noiseamp=0.0,noiseamp2=0.0,origAmp=0,ampOffset=0.0,switchLag=0.0,
		on=1,pos=0,minfreq=20,maxfreq=22000,minFil=100,recb,
		grainAmp=2.5,envb,inAmp2,thresh=0.135,fisConvAmp=0.75,
		fftb,fftb2,onsetsDiv=1,onsetProb=1.0,probs=#[1,1,1,1,1],convBuf,
		holdProb=0,holdDurMin=4,holdDurMax=5,envRel=1,silenceProb=0.0;

		var inAmp3,noise1,vcOrig,onsets,fft,fft2,inAmp,pitch,vc,in,
		imp, sound,freqs,generalenv,chain,grainOn,vcRec,envAtt,
		env,holdOnsets,holdOnsetsRel,holdDur,envSel,silence,convRecEnv,
		fisConv,convRecVc;

		generalenv = EnvGen.kr(Env.asr(0.5,1,6),on, doneAction: 2);
		vc = In.ar(inbus);
		vcOrig = vc;
		vc = Compander.ar(vc,vc,0.4,1,0.5) * 1.5;
		fft = FFT(LocalBuf(2048),vc);
		onsets = CoinGate.kr(onsetProb,
			PulseDivider.kr(Onsets.kr(fft,thresh.lag2(lagtime)),onsetsDiv));
		inAmp = Amplitude.ar(vc,0.01,0.2)*2;
		inAmp2 = ((Compander.ar(inAmp,inAmp,0.3,1,0.5) * 3) +
			ampOffset.lag2(lagtime)).clip(0.0,1.0);
		inAmp3 = (EnvFollow.ar(vc,0.99999) + ampOffset.lag2(lagtime)).clip(0.0,1.0);
		holdDur = TExpRand.kr(holdDurMin, holdDurMax, LocalIn.kr(1)+Impulse.kr(0.001));
		holdOnsets = SetResetFF.ar(K2A.ar(onsets),TDelay.ar(K2A.ar(onsets),holdDur));
		envAtt = TChoose.kr(holdOnsets,[0.01,0.01,0.01,0.5,1.5]);
		holdOnsetsRel = SetResetFF.ar(holdOnsets,TDelay.ar(holdOnsets,holdDur+envRel+envAtt));
		holdOnsets = holdOnsets * holdOnsetsRel;
		env = EnvGen.kr(Env.asr(envAtt,1,1),holdOnsets);
		envSel = TRand.ar(0,1,holdOnsets) < (holdProb.lag2(lagtime));
		silence = TRand.ar(0,1,holdOnsets) < (silenceProb.lag2(lagtime));

		onsets = Select.kr(envSel,[onsets,holdOnsets]);
		LocalOut.kr(onsets);
		inAmp2 = Select.kr(silence,[Select.kr(envSel, [inAmp2,env]),0]).lag2(0.01);
		inAmp3 = Select.kr(silence,[Select.kr(envSel, [inAmp3,env]),0]).lag2(0.01);
		amp = Select.kr(silence,[amp,0]).lag2(0.01);

		in = LocalIn.ar(1);
		in = Compander.ar(in,in,0.6,1,0.6) * 5;
		in = HPF.ar(in,30) * (inAmp > 0.005).lag2(0.2);
		in = in + (vc * 3);

		sound = SinOsc.ar(0,in*TExpRand.kr(
			rmin.lag2(lagtime),rmax.lag2(lagtime),onsets).lag2(2)) * 3;

		noise1 = FBSineN.ar(TExpRand.kr(2000,22050,onsets),
			TExpRand.kr(1.1,3,onsets),0.25,
			TExpRand.kr(1.15,5.5,onsets))*inAmp3;
		noise1 = WaveLoss.ar(noise1,TChoose.kr(onsets,[0,2,4,12]).lag2(1),40,2);
		noise1 = (LPCAnalyzer.ar(vc,noise1,512,16,windowtype:1) * 0.8)
		+ (noise1 * 0.3) + (vc * 0.05);
		noise1 = (noise1*4).softclip;
		noise1 = noise1 * Select.ar(TChoose.kr(onsets,[0,0,0,1]),[
			K2A.ar(1),Blip.ar(LFNoise1.ar(0.25).exprange(12,28)*
				LFNoise1.kr(3).range(0.75,1.5),3,0.25,0.75)]);

		sound = sound.softclip;
		sound = (LPCAnalyzer.ar(vc,sound,512,16,windowtype:1) * 0.8)
		+ (sound * 0.4);
		vcRec = Compander.ar(vc,vc,0.005,10,1);
		vcRec = Compander.ar(vcRec*2.25,vcRec*2.25,0.2,1,0.35)*1.25;
		RecordBuf.ar(vcRec,recb);

		fisConv = FBSineN.ar(TExpRand.kr(2000,22050,onsets).lag2(0.01),
			TExpRand.kr(1,4,onsets).lag2(3),0.25,
			TExpRand.kr(1.05,5.5,onsets).lag2(4)) * 0.5;
		fisConv = WaveLoss.ar(fisConv,TChoose.kr(onsets,[0,2,4,12,24]).lag2(0.5),40,2);
		convRecEnv = EnvGen.ar(Env.linen(0.01, 0.0264, 0.01), onsets);
		convRecVc = LPF.ar(HPF.ar(vcOrig,60),18000);
		convRecVc = (Compander.ar(convRecVc,convRecVc,0.2,1,0.3)*1.5)*convRecEnv;
		RecordBuf.ar(convRecVc,convBuf,loop:0,trigger:onsets);
		fisConv  = Convolution2L.ar(fisConv,convBuf,
			DelayC.kr(onsets,0.2,0.05),2048, crossfade: 1);
		fisConv = Compander.ar(fisConv,fisConv,0.3,1,0.5);
		fisConv = Compander.ar(fisConv,fisConv,0.3,1,0.4);
		fisConv = fisConv * TChoose.kr(onsets,[0.1,0.25,0.6,1.0,1.0]).lag2(TChoose.kr(onsets,[0.01,0.5,1,3]));

		imp = Impulse.ar(TExpRand.kr(12,40,onsets).lag2(2));
		sound = WaveLoss.ar(sound,TChoose.kr(onsets,[0,2,4,10]).lag2(1),40,2);
		sound = SelectX.ar(TChoose.kr(onsets,probs.lag2(lagtime).round).lag2(switchLag),
			[sound * noiseamp.lag2(2),
				((GrainBuf.ar(1,imp,TExpRand.kr(0.09,0.35,onsets),recb,
					1,TRand.kr(0.0,1.0,onsets).lag2(0.5),2,0,envb) * grainAmp.lag2(2) * 1.2)),
				noise1*noiseamp2.lag2(2),
				fisConv*fisConvAmp.lag2(2)]);

		sound = sound * TChoose.kr(onsets,[0.4,0.7,1.0,1.0]).lag2(TChoose.kr(onsets,[0.01,0.5,1,3]));
		sound = Compander.ar(sound,sound,0.3,1,0.5)*1.5;
		LocalOut.ar(FreeVerb.ar(sound,0.1,0.999)*(inAmp*4)+fbAmp.lag2(lagtime));
		sound = sound*amp.lag2(lagtime)*generalenv*(inAmp2+0.4);
		sound = Compander.ar(sound,sound,0.4,1,0.5);
		Out.ar(out,Pan2.ar(sound,pos) * 0.4);
	}).load;


	SynthDef("splitFold", {arg out=0,inbus=2,bufIn=0,busout,rmod=0.05,recBuf,window,
		lagtime=30,amp=0.1,rmin=3,rmax=16,analBuf,grainamp=0.1,
		on=1,minfreq=20,maxfreq=22000,minFil=100,onsProb=0.2,
		timeScaleMin=2,timeScaleMax=6,fbAmp=1,
		envAdd=1,thresh=0.16,envMinFr=0.1,envMaxFr=0.9,pan=0;
		var snd,trig,envs, onsets,fft,inAmp,pitch,vc,in,
		imp, sound,freqs,generalenv,onsetsEnv,sndRec;
		generalenv = EnvGen.kr(Env.asr(0.01,1,16),on, doneAction: 2);
		vc = In.ar(inbus);
		vc = Compander.ar(vc,vc,0.25,1,0.25) * 2.25;
		fft = FFT(LocalBuf(2048),vc);
		onsets = Onsets.kr(fft,thresh.lag2(lagtime));
		onsetsEnv = TRand.kr(0,1,onsets) < onsProb.lag2(lagtime);
		envs = EnvGen.kr(Env([0.001,0.3, 1,0.5,0.001],[0.8,0.8,1,1.75],[2,2,-2,-2]),onsetsEnv,
			timeScale:TExpRand.kr(timeScaleMin.lag2(lagtime),timeScaleMax.lag2(lagtime),onsetsEnv));
		inAmp = Amplitude.ar(vc,0.01,5);
		pitch = SpecCentroid.kr(fft).clip(100,4000).lag2(0.1);

		in = LocalIn.ar(1);
		in = Compander.ar(in,in,0.35,1,0.4) * 2;
		in = DynKlank.ar(
			`[[LFDNoise3.kr(0.005!4).exprange(minfreq,maxfreq),
				LFDNoise3.kr(0.005!4).exprange(minfreq,maxfreq)],
				[LFDNoise3.kr(0.01!4).exprange(0.05,0.3),
					LFDNoise3.kr(0.01!4).exprange(0.05,0.2)],
				[LFDNoise3.kr(0.01!4).exprange(0.25,3),
					LFDNoise3.kr(0.01!4).exprange(0.25,3)]
			], in
		) * 0.11 + in;

		in = (vc * 1.8) + (in * 1.9);

		rmin = (TRand.kr(0.5,1.25,PulseDivider.kr(onsets,5)).lag2(4) * rmin.lag2(lagtime));
		rmax = (TRand.kr(0.5,1.25,PulseDivider.kr(onsets,8)).lag2(4) * rmax.lag2(lagtime));
		sound = SinOsc.ar(0,in*LFNoise2.ar(rmod.lag2(0.1)).range(rmin.lag2(0.1),rmax.lag2(0.1)));
		sound = Mix.ar(sound);

		sound = LPCAnalyzer.ar(vc,sound,256,12,windowtype:0) * 0.6
		+ (sound * 0.5) * (inAmp * 5).clip(0.0,1.0);

		freqs = [ 77.387, 86.391, 118.59, 236.121, 459.544, 483.73, 806.099,
			1125.421, 2114.418, 4722.469, 7028.205, 8501.273, 11263.607 ];
		freqs = LFDNoise3.kr(0.07!13).exprange(freqs*0.75,freqs/0.75);
		sound = BHiPass4.ar(sound,[freqs[0],freqs[1],freqs[2],freqs[3],freqs[4],freqs[5],freqs[6],
			freqs[7],freqs[8],freqs[9],freqs[10],freqs[11]]);
		sound = BLowPass4.ar(sound,[freqs[1],freqs[2],freqs[3],freqs[4],freqs[5],freqs[6],
			freqs[7],freqs[8],freqs[9],freqs[10],freqs[11],freqs[12]]);

		sound = Compander.ar(sound,sound,0.3,1,0.1) * 2.6;
		sound = sound * envs;
		sound = sound * min(max(0, LFNoise1.kr(LFNoise1.kr(0.1!12).exprange(envMinFr.lag2(lagtime),
			envMaxFr.lag2(lagtime))) + envAdd.lag2(lagtime)),1.0).exprange;
		sound = BPeakEQ.ar(sound,pitch.lag2(0.4)*2,0.4,2);
		sound = BPeakEQ.ar(sound,pitch.lag2(0.4),0.3,3);
		sound = Compander.ar(sound,sound,0.1,1,0.5) * 2.1;
		sound = sound.softclip;
		sound = HPF.ar(sound,20);
		sndRec = (vc*1.3)+(Mix.ar(sound)*0.4);
		sndRec = Compander.ar(sndRec,sndRec,0.3,1,0.45)*2;

		RecordBuf.ar(sndRec,recBuf);
		trig = Impulse.ar(LFDNoise3.ar(LFDNoise3.ar(0.1).exprange(0.2,2.5)).exprange(13,40));
		snd = GrainBuf.ar(2,trig,LFNoise1.ar(0.1).exprange(0.05,0.22),
			recBuf,TRand.ar(0.95,1.05,K2A.ar(PulseDivider.kr(onsets,9))),
			LFNoise1.ar(0.1).range(0,1),
			2,TRand.ar(-1,1,trig),window,1024);
		snd = Compander.ar(snd,snd,0.4,1,0.25) * 2;
		snd = snd * EnvGen.kr(Env([0.001,1,1,0.001],[2,0.5,2.25],[2,1,-2]),PulseDivider.kr(onsets,6),
			timeScale:TExpRand.kr(0.75,4,PulseDivider.kr(onsets,6))) * (inAmp * 5).clip(0.1,1.0);

		sound = Mix.ar(sound);
		sound = sound.softclip;
		sound = FreeVerb.ar(sound,0.01,0.94);
		LocalOut.ar(sound + (Mix.ar(snd) * 0.1) * fbAmp.lag2(lagtime));

		Out.ar(out,(Pan2.ar(sound,pan.lag2(lagtime)) +
			(snd * grainamp.lag2(lagtime))) * amp.lag2(lagtime)*generalenv);
	}).load;



	SynthDef("unvoicedSec3", {arg out=0,inbus=2,bufIn=0,busout,
		rmod=0.05,lagtime=20,amp=0.1,rmin=3,rmax=16,
		noiseamp=0.0,noiseamp2=0.0,grainRates=#[1,1,1,1],envProbs=#[1,1,1,1,1,1],
		origAmp=0,ampOffset=0.0,switchLag=0.0,
		on=1,pos=0,minfreq=20,maxfreq=22000,minFil=100,recb,
		grainAmp=2.5,envb,inAmp2,thresh=0.135,convBuf,envAtt,
		fftb,fftb2,onsetsDiv=1,onsetProb=1.0,probs=#[1,1,1,1,1,1],lossGrAmp=1,
		holdProb=0,holdDurMin=4,holdDurMax=5,envRel=1,silenceProb=0.0;

		var inAmp3,noise1,vcOrig,onsets,fft,fft2,inAmp,pitch,vc,in,
		imp, sound,freqs,generalenv,chain,grainOn,vcRec,freq,
		env,holdOnsets,holdOnsetsRel,holdDur,envSel,silence,center,
		convRecEnv,convRecVc,demEnvs,lossGr;

		generalenv = EnvGen.kr(Env.asr(0.01,1,7),on, doneAction: 2);
		vc = In.ar(inbus);
		vcOrig = vc;
		vc = Compander.ar(vc,vc,0.4,1,0.5) * 1.5;
		fft = FFT(LocalBuf(2048),vc);
		center = SpecCentroid.kr(fft).clip(100,4000).lag2(0.1);

		onsets = CoinGate.kr(onsetProb,
			PulseDivider.kr(Onsets.kr(fft,thresh.lag2(lagtime)),onsetsDiv));
		inAmp = Amplitude.ar(vc,0.01,0.15)*2;
		inAmp2 = ((Compander.ar(inAmp,inAmp,0.3,1,0.5) * 3) +
			ampOffset.lag2(lagtime)).clip(0.0,1.0);
		inAmp3 = (EnvFollow.ar(vc,0.99999) + ampOffset.lag2(lagtime)).clip(0.0,1.0);

		holdDur = TExpRand.kr(holdDurMin.lag2(lagtime), holdDurMax.lag2(lagtime),
			LocalIn.kr(1)+Impulse.kr(0.001));
		holdOnsets = SetResetFF.ar(K2A.ar(onsets),TDelay.ar(K2A.ar(onsets),holdDur));
		envAtt = TChoose.kr(holdOnsets,[0.01,0.01,0.01,0.5,1.5])*envAtt;
		holdOnsetsRel = SetResetFF.ar(holdOnsets,TDelay.ar(holdOnsets,holdDur+envRel.lag2(lagtime)+envAtt));
		holdOnsets = holdOnsets * holdOnsetsRel;
		env = EnvGen.kr(Env.asr(envAtt,1,envRel.lag2(lagtime)),holdOnsets);
		envSel = TRand.ar(0,1,holdOnsets) < (holdProb.lag2(lagtime));
		silence = TRand.ar(0,1,holdOnsets) < (silenceProb.lag2(lagtime));
		onsets = Select.kr(envSel,[onsets,holdOnsets]);
		LocalOut.kr(onsets);
		inAmp2 = Select.kr(silence,[Select.kr(envSel, [inAmp2,env]),0]).lag2(0.01);
		inAmp3 = Select.kr(silence,[Select.kr(envSel, [inAmp3,env]),0]).lag2(0.01);

		amp = Select.kr(silence,[amp,0]).lag2(0.01);
		noise1 = FBSineN.ar(TExpRand.kr(500,22050,onsets),
			TExpRand.kr(1,4,onsets),0.25,
			TExpRand.kr(1.1,5.5,onsets))*inAmp3;
		noise1 = ((LPCAnalyzer.ar(vc,noise1,512,16,windowtype:1) * 0.8) )
		+ (noise1 * 0.3) + (vc * 0.05);
		noise1 = (noise1*4).softclip;
		noise1 = noise1 * Select.ar(TChoose.kr(onsets,[0,0,0,1]),[
			K2A.ar(1),Blip.ar(LFNoise1.ar(0.25).exprange(12,28)*
				LFNoise1.kr(3).range(0.75,1.5),3,0.25,0.75)]);

		vcRec = Compander.ar(vc,vc,0.005,10,1);
		vcRec = Compander.ar(vcRec*2.25,vcRec*2.25,0.2,1,0.35)*1.25;
		RecordBuf.ar(vcRec,recb);

		sound = FBSineN.ar(TExpRand.kr(2000,22050,onsets).lag2(0.01),
			TExpRand.kr(1,4,onsets).lag2(3),0.25,
			TExpRand.kr(1.05,5.5,onsets).lag2(4)) * 0.5;
		sound = WaveLoss.ar(sound,TChoose.kr(onsets,[0,2,4,12,24]).lag2(0.5),40,2);
		convRecEnv = EnvGen.ar(Env.linen(0.01, 0.0264, 0.01), onsets);
		convRecVc = LPF.ar(HPF.ar(vcOrig,60),18000);
		convRecVc = (Compander.ar(convRecVc,convRecVc,0.2,1,0.3)*1.5)*convRecEnv;
		RecordBuf.ar(convRecVc,convBuf,loop:0,trigger:onsets);
		sound  = Convolution2L.ar(sound,convBuf,
			DelayC.kr(onsets,0.2,0.05),2048, crossfade: 1);
		sound = Compander.ar(sound,sound,0.3,1,0.5);
		sound = Compander.ar(sound,sound,0.3,1,0.4);
		sound = sound * TChoose.kr(onsets,[0.1,0.25,0.6,1.0,1.0]).lag2(TChoose.kr(onsets,[0.1,0.5,1,3]));

		lossGr = DelayC.ar(vc,4,[0,0.8,2,3.1],[1.1,1,0.85,0.85]);
		lossGr = PitchShift.ar(lossGr,ExpRand(0.05!4,0.3),[1,0.85,1,0.8],0.1,0.1);
		lossGr = Mix.ar(WaveLoss.ar(lossGr, Rand(72,86!4), 90,2));
		lossGr = Compander.ar(lossGr,lossGr, 0.1,1,0.2)*3;

		imp = Impulse.ar(TExpRand.kr(10,40,onsets).lag2(2));
		sound = SelectX.ar(TChoose.kr(onsets,probs.lag2(lagtime).round).lag2(switchLag.lag2(lagtime)),
			[sound * noiseamp.lag2(lagtime),
				((GrainBuf.ar(1,imp,TExpRand.ar(0.09,0.35,K2A.ar(onsets)),recb,
					TChoose.ar(K2A.ar(onsets),K2A.ar(grainRates)),
					TRand.ar(0.1,0.9,K2A.ar(onsets)).lag2(0.5),
					2,0,envb) * grainAmp.lag2(lagtime))), noise1*noiseamp2.lag2(lagtime),
				lossGr * lossGrAmp.lag2(lagtime)]);

		sound = sound * SelectX.kr(TChoose.kr(onsets,envProbs.lag2(1)).round.lag2(0.05),
			[1,EnvGen.kr(Env([0.001,1,1,0.001],[0.4,0.07,0.05]),onsets,timeScale:TRand.kr(0.5,2,onsets)),
				EnvGen.kr(Env([0.001,1,1,0.001],[0.01,0.65,0.02]),onsets,timeScale:TRand.kr(0.5,2,onsets)),
				EnvGen.kr(Env([0.001,1,1,0.001],[0.01,0.1,0.6]),onsets,timeScale:TRand.kr(0.5,2,onsets))]);
		sound = LeakDC.ar(sound);
		sound = Compander.ar(sound,sound,0.3,1,0.5)*1.5;
		sound = sound*amp.lag2(lagtime)*generalenv*((inAmp2+0.35).clip(0.35,1.0));
		sound = Compander.ar(sound,sound,0.4,1,0.5);
		Out.ar(out,Pan2.ar(sound,pos) * 0.75);
	}).load;

	SynthDef("pitchFreeze5", {arg out=0,inBus,holdDurMin=2,holdDurMax=10,
		lagtime=1,trigsProb=0.06,
		amp=0.23,recB,on=1;
		var vc, sndRec,gr,pitch,pitches,hasFreq,sins,inp,del,trigs,n=5,
		holdOnsets,holdDur,holdOnsetsRel,envAtt,trigsOrig,grTrig,generalenv;
		vc = In.ar(inBus);
		# pitch, hasFreq = Pitch.kr(vc*1.5,440,60,4000,peakThreshold:0.5);
		sndRec = Compander.ar(vc,vc,0.6,1,0.6)*1.5;

		generalenv = EnvGen.kr(Env.asr(0.01,1,7),on, doneAction: 2);

		trigsOrig = Onsets.kr(FFT(LocalBuf(2048),vc), 0.13) * hasFreq;
		trigs = TRand.kr(0.0,1.0,trigsOrig!n) < trigsProb.lag2(lagtime);
		trigs = trigs * hasFreq;

		holdDur = TExpRand.kr(holdDurMin.lag2(lagtime), holdDurMax.lag2(lagtime),
			LocalIn.kr(n)+Impulse.kr(0.001));
		holdOnsets = SetResetFF.ar(K2A.ar(trigs),TDelay.ar(K2A.ar(trigs),holdDur));
		//envAtt = TChoose.kr(holdOnsets,[1,2,2.5,3.6,4,6,7]);
		envAtt = TExpRand.kr(6,16,holdOnsets);
		holdOnsetsRel = SetResetFF.ar(holdOnsets,TDelay.ar(holdOnsets,holdDur+envAtt));
		holdOnsets = holdOnsets * holdOnsetsRel;
		LocalOut.kr(holdOnsets);

		pitches = K2A.ar(Latch.kr(pitch!n, holdOnsets).clip(10,1000));

		sins = PhinGen3.ar(1, // 1 reader
			K2A.ar(300), // stiffness
			K2A.ar(200), K2A.ar(150), // dampings
			K2A.ar(50000), // mass
			K2A.ar(4), // control points
			K2A.ar(1), // distance
			LFDNoise1.ar(0.1!n).exprange(1,30), // updateFreq
			// control/excite the model
			[
				// the reader (frequency)
				pitches*LFNoise1.ar(0.03).range(0.9,1.1),

				SinOsc.ar(ExpRand(0.02!n,0.25)),
				SinOsc.ar(ExpRand(0.02!n,0.25)),
				SinOsc.ar(ExpRand(0.02!n,0.25)),
				SinOsc.ar(ExpRand(0.02!n,0.25))
		]) * 0.4;
		sins = sins * EnvGen.kr(Env([0,0,1,1],[4,2,500]),1);
		sins = sins *
		TChoose.kr(holdOnsets, [
			EnvGen.kr(Env([0,1,1,0],[0.5,0.09,0.41],[3,-2,-2]),holdOnsets,timeScale:holdDur),
			EnvGen.kr(Env([0,1,1,0],[0.4,0.1,0.5],[3,-2,-2]),holdOnsets,timeScale:holdDur),
			EnvGen.kr(Env([0,1,1,0],[0.45,0.05,0.5],[3,-2,-2]),holdOnsets,timeScale:holdDur),
			EnvGen.kr(Env([0,1,1,0],[0.3,0.1,0.6],[3,-2,-2]),holdOnsets,timeScale:holdDur)
		]);

		RecordBuf.ar(sndRec+(Mix.ar(sins)*0.0),recB);
		gr = TGrains.ar(2,Impulse.ar(LFNoise1.ar(0.1).exprange(10,60)),recB,1,
			LFNoise1.ar(0.07).range(0.5,3.5),0.1) * 2;
		grTrig = (TRand.kr(0.0,1.0,trigsOrig) < 0.1);
		gr = gr * EnvGen.kr(Env([0,1,1,0],[1.5,1,1.5],[2,-2,-2]),grTrig,timeScale:TRand.kr(0.6,1.25,grTrig));

		Out.ar(out, ((gr*0) + FreeVerb.ar(Splay.ar(sins),0.1,0.8)) * amp.lag2(lagtime) * generalenv);
	}).load(s);


	SynthDef("unvoiced_section4", {arg out = 0, amp = 0.5, on = 1, inbus, recb;
		var vc, sndRec,gr,pitch,pitches,hasFreq,sins,inp,del,trigs,n=8,generalenv,envProbs=#[1,2,3,4],
		holdOnsets,holdDur,holdOnsetsRel,holdDurMin,holdDurMax,envAtt,trigsOrig,grTrig,vcOnsets;
		vc = In.ar(inbus);
		vc = Compander.ar(vc,vc,0.2,1,0.3)*2;
		vc = FreeVerb.ar(vc,0.15,0.8);
		# pitch, hasFreq = Pitch.kr(vc*1.5,440,60,4000,peakThreshold:0.6);
		sndRec = Compander.ar(vc,vc,0.4,1,0.5)*1.6;
		holdDurMin = 2.25;
		holdDurMax = EnvGen.kr(Env([3.5,4,10,10,5],[15,30,30,20]),1);

		generalenv = EnvGen.kr(Env.asr(0.01,1,7),on, doneAction: 2);
		vcOnsets = Onsets.kr(FFT(LocalBuf(2048),vc), 0.08);
		trigsOrig = vcOnsets + Dust.kr(0.8) * hasFreq;
		trigs = TRand.kr(0.0,1.0,trigsOrig!n) < 0.8;
		trigs = trigs + (Dust.kr(0.25!n) * hasFreq);

		holdDur = TExpRand.kr(holdDurMin, holdDurMax, LocalIn.kr(n)+Impulse.kr(0.001));
		holdOnsets = SetResetFF.ar(K2A.ar(trigs),TDelay.ar(K2A.ar(trigs),holdDur));
		envAtt = TExpRand.kr(EnvGen.kr(Env([9,7,3,1,2],[15,25,30,25]),1),12,holdOnsets);
		holdOnsetsRel = SetResetFF.ar(holdOnsets,TDelay.ar(holdOnsets,holdDur+envAtt));
		holdOnsets = holdOnsets * holdOnsetsRel;
		LocalOut.kr(holdOnsets);

		pitches = K2A.ar(Latch.kr(pitch!n, holdOnsets).clip(10,3000));
		holdOnsets = holdOnsets * (pitches > 240);

		sins = PhinGen3.ar(1, // 1 reader
			K2A.ar(300), // stiffness
			K2A.ar(200), K2A.ar(150), // dampings
			K2A.ar(50000), // mass
			K2A.ar(4), // control points
			K2A.ar(1), // distance
			LFDNoise1.ar(0.1!n).exprange(1,30), // updateFreq
			// control/excite the model
			[
				// the reader (frequency)
				pitches*LFNoise1.ar(0.04).range(0.97,1.03),
				//*Rand(0.95!n,1.05),//
				SinOsc.ar(ExpRand(0.02!n,0.25)),
				SinOsc.ar(ExpRand(0.02!n,0.25)),
				SinOsc.ar(ExpRand(0.02!n,0.25)),
				SinOsc.ar(ExpRand(0.02!n,0.25))
		]) * 0.4;
		sins = sins * EnvGen.kr(Env([0,0,1,1],[4,2,500]),1);
		sins = sins *
		SelectX.kr(TChoose.kr(holdOnsets, envProbs), [
			EnvGen.kr(Env([0,1,1,0],[0.02,0.25,0.73],[3,-2,-2]),holdOnsets,timeScale:holdDur),
			EnvGen.kr(Env([0,1,1,0],[0.25,0.05,0.7],[3,-2,-2]),holdOnsets,timeScale:holdDur),
			EnvGen.kr(Env([0,1,1,0],[0.3,0.1,0.6],[3,-2,-2]),holdOnsets,timeScale:holdDur),
			EnvGen.kr(Env([0,1,1,0],[0.4,0.1,0.5],[3,-2,-2]),holdOnsets,timeScale:holdDur),
			EnvGen.kr(Env([0,1,1,0],[0.45,0.05,0.5],[3,-2,-2]),holdOnsets,timeScale:holdDur)
		]);

		RecordBuf.ar((sndRec*1.2)+(Mix.ar(sins)*0.2),recb);
		gr = TGrains.ar(2,Impulse.ar(LFNoise1.ar(0.1).exprange(10,60)),recb,1,
			LFNoise1.ar(0.07).range(0.25,3.5),0.1) * 2;
		grTrig = (TRand.kr(0.0,1.0,trigsOrig) < 0.6) + Dust.kr(0.2);
		gr = gr * EnvGen.kr(Env([0,1,1,0],[1.5,1,1.5],[2,-2,-2]),grTrig,timeScale:TRand.kr(0.6,1.25,grTrig));

		Out.ar(out, (gr*2) + FreeVerb.ar(Splay.ar(sins),0.12,0.85) * amp.lag2(3) * 0.8 * generalenv);
	}).load;


	// SECTIONS

	// 1
	~section1 = Task({
		~recb1.zero;
		~recb2.zero;
		~recb3.zero;
		~recb4.zero;
		~analBuf1.zero;
		~analBuf2.zero;
		~analBuf3.zero;
		a = Synth.after(~inputPlay,"unvoicedSec1", [\out, ~output, \minFil, 30, \rmod, 0.05,
			\rmin, 0.75, \rmax, 10, \amp, 0.42, \pos, -1, \inbus, ~voiceBus.index,
			\thresh, 0.23,\fftb,~b1.bufnum,\fftb2,~b2.bufnum,\envAdd, 1.7,
			\recb,~recb1,\delay,0,\grainRates,#[1,1,2,1],\origAmp,0,
			\noiseamp,0.8,\grainAmp,0.9,\envb,~wind.bufnum,\probs,#[1,1,1,1,1],\convBuf,~analBuf1.bufnum]);
		b =  Synth.after(~inputPlay,"unvoicedSec1", [\out, ~output, \minFil, 30, \rmod, 0.05,
			\rmin, 0.75, \rmax, 10, \amp, 0.42, \pos, 1, \inbus, ~voiceBus.index, \thresh, 0.23,
			\fftb,~b3.bufnum,\fftb2,~b4.bufnum,\envAdd, 1.7,\recb,~recb2,\delay,0,
			\grainRates,#[1,1,2,1],\origAmp,0,\noiseamp,0.8,\grainAmp,0.9,
			\envb,~wind.bufnum,\probs,#[1,1,1,1,1],\convBuf,~analBuf2.bufnum]);

		39.wait;
		a.set(\probs,#[0,0,1,1,1],\grainRates,#[0.5,1,1,2],\noiseamp,0.21,\onsetProb,0.6,
			\amp,0.52,\holdProb, 0.12,\lagtime, 16,\silenceProb, 0.1,\holdDurMin, 1, \holdDurMax, 4);
		b.set(\probs,#[0,0,1,1,1],\grainRates,#[0.5,1,1,2],\noiseamp,0.21,\onsetProb,0.6,
			\amp,0.52,\holdProb, 0.12,\lagtime, 18,\silenceProb, 0.1,\holdDurMin, 1, \holdDurMax, 4);


		39.wait; //78
		a.set(\lagtime,0.1,\probs,#[0,0,3,3,3],\onsetProb,1,\grainRates,#[1,1,1,1],
			\noiseamp,0.22,\rmin, 1, \rmax, 10, \noiseamp2,0.23,\envAdd, 1,\delay,0,
			\thresh,0.18,\amp, 0.57,\holdProb, 0.6, \holdDurMin, 2, \holdDurMax, 6,
			\envRel, 0.1, \silenceProb, 0.24);
		b.set(\lagtime,0.1,\probs,#[0,0,3,3,3],\onsetProb,1,\grainRates,#[1,1,1,1],\noiseamp,0.22,
			\rmin, 1, \rmax, 10, \noiseamp2,0.23, \envAdd, 1,\thresh,0.18,\delay,0,
			\amp, 0.57,\holdProb, 0.6,  \holdDurMin, 2, \holdDurMax, 6, \envRel, 0.1, \silenceProb, 0.24);

		25.5.wait; // 103.5
		w = Synth.after(~inputPlay,"unvoicedSec1",[\out, ~output, \minFil, 30, \rmod, 0.05, \rmin, 1.2,
			\rmax, 16, \amp, 0.55, \pos, 0, \inbus, ~voiceBus.index, \thresh, 0.2,\fftb,~b5.bufnum,
			\fftb2,~b6.bufnum,\envAdd, 1.5,\recb,~recb3.bufnum,\delay,0,\grainRates,#[1,1,1,1],
			\origAmp,0,\onsetProb,0.9,\probs,#[1,2,2,2,3],\noiseamp,0.25,\noiseamp2,0.25,
			\grainAmp,0.8,\convBuf,~analBuf3.bufnum,\envb,~wind.bufnum,
			\holdProb, 0.6, \silenceProb, 0.12, \holdDurMin, 3, \holdDurMax, 8, \envRel, 3, \lagtime, 15]);
		a.set(\lagtime,1.1,\silenceProb, 0.15,\onsetProb,0.8,  \probs,#[1,1,1,2,3], \amp, 0.55);
		b.set(\lagtime,1.1,\silenceProb, 0.15, \onsetProb,0.8, \probs,#[1,1,1,2,3], \amp, 0.55);

		13.5.wait; // 117
		a.set(\lagtime, 20, \ampOffset,0.7,\onsetProb,0.8,\amp, 0.62,\switchLag,3,\rmax,4,
			\probs,#[1,1,2,2,3],\noiseamp2,0.26,\holdProb, 0.95, \holdDurMin, 3, \holdDurMax, 9,
			\envRel, 0.1, \silenceProb, 0.0,\grainAmp,0.8);
		b.set(\lagtime,20, \ampOffset,0.7,\onsetProb,0.8,\amp, 0.62,\switchLag,4,\rmax,4,
			\probs,#[1,1,2,2,3],\noiseamp2,0.26,\holdProb, 0.95, \holdDurMin, 3, \holdDurMax, 9,
			\envRel, 0.1, \silenceProb, 0.0,\grainAmp,0.8);
		w.set(\lagtime,20, \ampOffset,0.7,\amp, 0.6,\onsetProb,0.8,\switchLag,5,\rmax,6,
			\probs,#[1,1,2,2,3],\noiseamp2,0.26,\envRel, 4,\holdProb, 0.95, \holdDurMin, 3,
			\holdDurMax, 9, \silenceProb, 0.0);

		17.wait;
		// play file
		{ Out.ar(~output, PlayBuf.ar(2, ~trans.bufnum,doneAction: 2) * 0.7) }.play;

		2.wait; // 136
		a.set(\lagtime, 6, \fbAmp,0.7, \ampOffset,1,\onsetProb,0.65,\amp, 0.67,\switchLag,4,
			\rmax,4,\probs,#[0,1,2,2,2],\noiseamp2,0.3,\holdProb, 0.9, \holdDurMin, 7, \holdDurMax, 12,
			\envRel, 0.25, \silenceProb, 0.0,\grainAmp,0.8);
		b.set(\lagtime,7, \fbAmp,0.7, \ampOffset,1,\onsetProb,0.65,\amp, 0.67,\switchLag,3,\rmax,4,
			\probs,#[0,1,2,2,2],\noiseamp2,0.3,\holdProb, 0.9, \holdDurMin, 6, \holdDurMax, 12,
			\envRel, 0.25, \silenceProb, 0.0,\grainAmp,0.8);
		w.set(\lagtime,8, \fbAmp,0.7, \ampOffset,1,\amp, 0.62,\switchLag,5,\rmax,6,
			\probs,#[0,2,2,2,3],\noiseamp2,0.3,\envRel, 0.25,\holdProb, 0.9, \holdDurMin, 6,
			\holdDurMax, 12, \silenceProb, 0.0);


		18.wait;
		a.set(\on, 0);
		b.set(\on, 0);
		w.set(\on, 0);
	});

	// 2
	~section2 = Task({
		x = Synth.after(~inputPlay,"splitFold", [\out, ~output, \minFil, 70, \rmod, 0.1,
			\rmin, 0.07, \rmax, 15, \amp, 0.15,  \inbus, ~voiceBus.index, \envAdd,0.4,
			\analBuf, ~analBuf1.bufnum,\recBuf,~recS21.bufnum, \window, ~winBuf.bufnum,\pan,-1]);
		y = Synth.after(~inputPlay,"splitFold", [\out, ~output, \minFil, 70, \rmod, 0.1,
			\rmin, 0.07, \rmax, 15, \amp, 0.15,  \inbus, ~voiceBus.index, \envAdd,0.4,
			\analBuf, ~analBuf2.bufnum,\recBuf,~recS22.bufnum, \window, ~winBuf.bufnum,\pan,1]);

		15.wait; // 174
		x.set(\rmax, 5.5, \lagtime, 14, \pan, -1, \amp, 0.09, \fbAmp, 0.7);
		y.set(\rmax, 5.5, \lagtime, 14, \pan, 1, \amp, 0.09, \fbAmp, 0.7);

		18.75.wait; // 192.75
		x.set(\rmax, 2.0, \lagtime, 40, \pan, -1, \grainamp, 0.72, \amp, 0.15, \fbAmp, 1);
		y.set(\rmax, 2.0, \lagtime, 40, \pan, 1, \grainamp, 0.72, \amp, 0.15, \fbAmp, 1);

		52.5.wait; // 245.25
		x.set(\rmax, 2.0, \timeScaleMin, 1.3, \timeScaleMax, 4.0, \onsProb, 0.15,
			\lagtime, 16, \fbAmp, 0.65, \pan, -1, \amp, 0.12);
		y.set(\rmax, 2.0, \timeScaleMin, 1.3, \timeScaleMax, 4.0, \onsProb, 0.15, \lagtime, 16,
			\fbAmp, 0.65, \pan, 1, \amp, 0.12);

		19.75.wait; //265
		x.set(\on, 0, \grainamp, 0.12, \lagtime, 7);
		y.set(\on, 0, \grainamp, 0.12, \lagtime, 7);
	});

	~section3 = Task({
		~recb1.zero;
		~recb2.zero;
		~recb3.zero;
		~recb4.zero;
		~analBuf1.zero;
		~analBuf2.zero;
		~analBuf3.zero;
		~recS31.zero;
		a = Synth.after(~inputPlay,"unvoicedSec3", [\out, ~output, \ampOffset,0.5,\onsetProb,0.7,\amp, 0.35,
			\switchLag,0.1,\rmod, 0.08, \rmin, 1, \rmax,15,\probs,#[1,2,1,0,0,3],\holdProb, 0.15,
			\holdDurMin, 2,\holdDurMax, 9, \envRel, 2, \silenceProb, 0.1,\grainAmp,0.75,
			\lossGrAmp, 1.7, \pos, -1, \inbus, ~voiceBus.index,
			\thresh, 0.2,\fftb,~b1.bufnum,\fftb2,~b2.bufnum,\envAdd, 1.75,
			\recb,~recb1,\convBuf,~analBuf1.bufnum,
			\delay,0,\grainRates,#[1,1,1,1],\origAmp,0,\noiseamp,1.2,\noiseamp2,0.22,
			\envb,~wind.bufnum,\envProbs,#[0,0,1,0,1,2]]);
		b = Synth.after(~inputPlay,"unvoicedSec3", [\out, ~output, \ampOffset,0.5,\onsetProb,0.7,\amp, 0.35,
			\switchLag,0.1,\rmod, 0.08,\rmin, 1, \rmax,15 ,\probs,#[1,2,1,0,0,3],\holdProb, 0.15,
			\holdDurMin, 2, \holdDurMax, 9,\envRel, 2, \silenceProb, 0.1,\grainAmp,0.75,
			\lossGrAmp, 1.7, \pos, 1, \inbus, ~voiceBus.index,
			\thresh, 0.2,\fftb,~b3.bufnum,\fftb2,~b4.bufnum,\envAdd, 1.75,
			\recb,~recb2,\convBuf,~analBuf2.bufnum,\delay,0,\grainRates,#[1,1,1,1],
			\origAmp,0,\noiseamp,1.2,\noiseamp2,0.22,\envb,~wind.bufnum,
			\envProbs,#[0,0,1,0,1,2]]);
		p = Synth.after(~inputPlay,"pitchFreeze5", [\out, ~output, \inBus, ~voiceBus.index,
			\recB, ~recS31.bufnum, \trigsProb,0.04, \amp, 0.0]);
		c = Synth.after(~inputPlay,"unvoicedSec3", [\out, ~output, \ampOffset,0.5,\onsetProb,0.7,\amp, 0.35,
			\switchLag,0.1,\rmod, 0.08, \rmin, 1, \rmax,15 ,\probs,#[1,2,1,0,0,3],
			\holdProb, 0.15, \holdDurMin, 2, \holdDurMax, 9, \envRel, 2, \silenceProb, 0.1,
			\grainAmp,0.75, \lossGrAmp, 1.7, \pos, 0, \inbus, ~voiceBus.index,
			\thresh, 0.2,\fftb,~b5.bufnum,\fftb2,~b6.bufnum,\envAdd, 1.75,
			\recb,~recb3,\convBuf,~analBuf3.bufnum, \delay,0,\grainRates,#[1,1,1,1],
			\origAmp,0,\noiseamp,1.2,\noiseamp2,0.22,\envb,~wind.bufnum,
			\envProbs,#[0,0,1,0,1,2]]);

		20.wait; // 298
		a.set(\lagtime, 30, \switchLag, 0.1, \rmod, 0.08, \rmin, 1, \rmax, 6,
			\amp, 0.35, \grainRates,#[1,1,1,1],\probs,#[1,2,1,0,0,3],\holdProb, 0.5,
			\holdDurMin, 2, \holdDurMax, 10, \envRel, 2.25, \silenceProb, 0.1,\envProbs,#[0,0,1,0,0,2]);
		b.set(\lagtime, 36, \switchLag, 0.1, \rmod, 0.08, \rmin, 1, \rmax, 6,
			\amp, 0.35, \grainRates,#[1,1,1,1],\probs,#[1,2,1,0,0,3],\holdProb, 0.5,
			\holdDurMin, 2, \holdDurMax, 10, \envRel, 2.25, \silenceProb, 0.1, \envProbs,#[0,0,1,0,0,2]);
		c.set(\lagtime, 32, \switchLag, 0.1,  \rmod, 0.08, \rmin, 1, \rmax, 6,
			\amp, 0.35, \grainRates,#[1,1,1,1],\probs,#[1,2,1,0,0,3],\holdProb, 0.5,
			\holdDurMin, 2, \holdDurMax, 10, \envRel, 2.25, \silenceProb, 0.1, \envProbs,#[0,0,1,0,0,2]);
		p.set(\lagtime, 25, \amp, 0.28);

		14.wait; // 312
		a.set(\lagtime, 24, \switchLag, 0.2, \probs,#[1,1,0,0,0,3], \grainAmp,0.68, \envProbs,#[0,0,1,0,0,2]);
		b.set(\lagtime, 24, \switchLag, 0.2, \probs,#[1,1,0,0,0,3], \grainAmp,0.68, \envProbs,#[0,0,1,0,0,2]);
		p.set(\lagtime, 20, \amp, 0.37, \trigsProb, 0.07);

		16.wait; // 328
		a.set(\amp, 0.32, \lagtime, 24, \probs,#[1,1,0,0,0,0], \holdProb, 0.75, \grainAmp,0.6, \silenceProb, 0.25,
			\holdDurMin, 4, \holdDurMax, 10, \noiseamp, 1.1, \envRel, 2.5, \envAtt, 3,
			\envProbs,#[0,0,0,0,0,0]);
		b.set(\amp, 0.32, \lagtime, 24, \probs,#[1,1,0,0,0,0], \holdProb, 0.75, \grainAmp,0.6, \silenceProb, 0.25,
			\holdDurMin, 4, \holdDurMax, 10, \noiseamp, 1.1, \envRel, 2.5, \envAtt, 3,
			\envProbs,#[0,0,0,0,0,0]);
		c.set(\lagtime, 4, \probs,#[1,2,1,0,0,0], \silenceProb, 0.15);

		8.wait; // 336
		p.set(\lagtime, 8, \amp, 0.5, \trigsProb, 0.12, \holdDurMin, 4);
		a.set(\amp, 0.28, \lagtime, 6, \silenceProb,0.36, \probs,#[1,1,0,0,0,1], \envProbs,#[0,0,0,0,0,0]);
		b.set(\amp, 0.28, \lagtime, 6, \silenceProb, 0.36, \probs,#[1,1,0,0,0,1], \envProbs,#[0,0,0,0,0,0]);
		c.set(\amp, 0.28, \lagtime, 7, \holdProb, 0.75, \grainAmp,0.6, \silenceProb, 0.3,
			\holdDurMin, 4, \holdDurMax, 10, \noiseamp, 1.1, \envRel, 2.5, \envAtt, 3,
			\silenceProb, 0.36, \probs,#[1,1,0,0,0,1], \envProbs,#[0,0,0,0,0,0]);

		24.wait; // 360
		a.set(\on, 0);
		b.set(\on, 0);
		c.set(\on, 0);
		p.set(\on, 0);

	});

	~section4 = Task({
		~recS4.zero;
		x = Synth.after(~inputPlay,"unvoiced_section4", [\amp, 0.6, \inbus, ~voiceBus.index,
			\recb, ~recS4, \envProbs, #[1,2,3,4], \out, ~output]);
		80.wait; // 442
		x.set(\envProbs, #[0,0,0,2]);
		28.wait; // ending
		x.set(\on, 0);


	});

	// PIECE
	~unvoiced = Task({
		"section 1".postln;
		~section1.start;
		159.wait;

		"section 2".postln;
		~section2.start;
		119.wait;

		"section 3".postln;
		~section3.start;
		84.wait;

		"section 4".postln;
		~section4.start;
	});
});
)




// testing with recording
//~allsec = Buffer.read(s, "/home/luc/Work/voice_elec/sec3.wav");


// RUN FIRST
(
~inputPlay = {
	var amp,snd;
	// testing with recording
	//snd = PlayBuf.ar(1,~allsec.bufnum,loop: 0,startPos:44100*0);
	snd = SoundIn.ar(~input);
	amp = 0.9; // level of direct amplification, adjust!
	snd = Compander.ar(snd,snd,0.2,1,0.3)*2.5;
	snd = Compander.ar(snd,snd,0.2,1,0.3);
	snd = FreeVerb.ar(snd,0.06,0.77);
	Out.ar(~voiceBus.index,snd*2);
	Out.ar(~outputVc,snd.dup*amp);
}.play;
)


// START LIVE-ELECTRONICS
~unvoiced.start;


// RECORDING

//s.recChannels = 4;
//s.prepareForRecord("/home/luc/Work/voice_elec/recberlin1.wav");
//s.record;
//s.stopRecording;


// for rehearsals

//v = s.volume;
//v.gui;

~section1.start;

~section2.start;

~section3.start;

~section4.start;

s.scope

s.reboot;

s.meter