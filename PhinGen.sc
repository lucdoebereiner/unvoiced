
PhinGen3 : MultiOutUGen {
	*ar {
		arg nReaders = 1, stiffness = 200, viscosity = 200, d2 = 20, mass = 50000, n = 10, distance = 10, updateFreq = 300, array;
		// distance can never grow above the initial value
		^this.multiNewList(['audio', nReaders, stiffness, viscosity, d2, mass, n, distance, updateFreq] ++ array)
	}

	init {arg ... theInputs;
		inputs = theInputs;
		^this.initOutputs(inputs[0], rate);
	}
}





