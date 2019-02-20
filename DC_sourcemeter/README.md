# DC_sourcemeter
A tool for making DC measurements with a sourcemeter

## Usage
```bash
$ ./run_me.py --help
usage: run_me.py [-h] -a SM_ADDRESS [-v] -s SET_POINT -t DURATION -c
                 COMPLIANCE [-n NPLC] [-f] [-w]

Make DC sourcemeter measurements under constant current or voltage and dump
csv to stdout: time, current, voltage, status

optional arguments:
  -h, --help            show this help message and exit
  -a SM_ADDRESS, --sm-address SM_ADDRESS
                        Sourcemeter connection address, eg. "ASRL0::INSTR"
  -v, --source-voltage  Source voltage (if this argument is absent current
                        will be sourced)
  -s SET_POINT, --set-point SET_POINT
                        Source value in amps or volts
  -t DURATION, --duration DURATION
                        Number of seconds to measure for
  -c COMPLIANCE, --compliance COMPLIANCE
                        Compliance value in amps or volts
  -n NPLC, --nplc NPLC  NPLC value
  -f, --front           Use the front terminals
  -w, --four-wire       Use four wire mode
```
## Example usage
```bash
$ ./run_me.py --sm-address ASRL/dev/ttyUSB0::INSTR --duration 10 --set-point 0.001 --compliance 2 --nplc 10 > data.csv
```
- for serial communication with a sourcemeter attached to /dev/ttyUSB0
- sourcing 1mA
- for 10 seconds
- with a voltage compliance setting of 2V
- with NPLC=10 integration
- and storing the data generated in a file called data.csv

## Example output
```bash
$ cat data.csv 
#Time [s],Current [A],Voltage [V],Status Bits
0.623,-0.000000,2.000040,39944
1.068,-0.000000,2.000039,39944
1.513,-0.000000,2.000039,39944
1.957,-0.000000,2.000039,39944
2.400,-0.000000,2.000040,39944
2.844,-0.000000,2.000040,39944
3.288,-0.000000,2.000040,39944
3.732,-0.000000,2.000041,39944
4.176,-0.000000,2.000041,39944
4.619,-0.000000,2.000041,39944
5.062,-0.000000,2.000041,39944
5.506,-0.000000,2.000040,39944
5.950,-0.000000,2.000041,39944
6.395,-0.000000,2.000040,39944
6.838,-0.000000,2.000040,39944
7.281,-0.000000,2.000040,39944
7.725,-0.000000,2.000040,39944
8.168,-0.000000,2.000041,39944
8.611,-0.000000,2.000041,39944
9.056,-0.000000,2.000041,39944
9.499,-0.000000,2.000041,39944
9.943,-0.000000,2.000041,39944
```
Then you could read that data in with the pandas python module like so:
```python
import pandas
the_data = pandas.read_csv('data.csv')
```
