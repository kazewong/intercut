The use this program, the following package needed to be installed:

	numpy
	scipy
	astropy
	mpmath

The following files are needed as the input of the program:

- A noise file,specifying the noise correspoding to a specific wavelength. If radial dependent noise is needed, please set the first value of flags to be 1

The format of the noise file should follow:

           (without radial dependent noise)

 	wavelength			noise

            (with radial dependent noise)

				radius of the object

	wavelength			noise


- A parameter file, specifying the name and the wavelength of the lines. Please ensure if the name of lines is same as the name in the catalogue. Note that the first row is reserved for the main line, please enter the parameter of your main line there, the order of the rest does not matter.

- A condition file, specifying all the conditions that the lines will need to pass.Predefined variable can be used to simplified the input , those variables can be found in the section predefined variable in the program. 

The following input will be optional:

- Photometric cut for the whole catalogue,predefined variables can also be used in this section

- Recalibration of lines using different luminosity function.(This part require users to write their own luminosity function, please look into the LuminosityFunction.py , the function CalibrationHa and CalibrationO3b is the predefined luminosity function, customized luminosity function could be defined in a similar way.)

- If fluctuation is needed in flux, please set the second values of the flags to be 1
 
The way to use the control script:

each row of the control script should have two component, the first component specify what kind of argument it is, and the second argument is the argument itself,please seperate the two components by ;

flag: on/off switch for random fluctuation and radial dependent noise
catalog: input catalogue
noise: input noise
lines: lines paramter file
output: desired output text file
loopNumber: simluation number, the second line identification part can be looped for multiple times, which the flux should be affected by random fluuctuation
binNumber: output histogram bin number
lowerLim: output histogram lower x limit
upperLim: output histogram upper x limit
photocut: photometric cut for the whole survey
condition: condition file
noisepara: this is needed only when radial noise is needed, the first 3 arguments are the range and bin of wavelength, and the last 3 arguments are range and bin of radius of the object.



