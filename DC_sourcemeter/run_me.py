#!/usr/bin/env python3

# written by grey@christoforo.net
# on 20 Feb 2019

# Example usage:
# ./run_me.py --sm-address ASRL/dev/ttyUSB0::INSTR --duration 10 --set-point 0.001 --compliance 2 --nplc 10 > data.csv
# for serial communication with a sourcemeter attached to /dev/ttyUSB0
# sourcing 1mA
# for 10 seconds
# with a voltage compliance setting of 2V
# with NPLC=10 integration
# and storing the data generated in a file called data.csv

# data generated here can be read in nicely with pandas:
# pandas.read_csv('data.csv')

from k2400 import k2400
import argparse

parser = argparse.ArgumentParser(description='Make DC sourcemeter measurements under constant current or voltage and dump csv to stdout: time, current, voltage, status')
parser.add_argument('-a', '--sm-address', type=str, required=True, help='Sourcemeter connection address, eg. "ASRL0::INSTR"')
parser.add_argument('-v', '--source-voltage', default=False, action='store_true', help="Source voltage (if this argument is absent current will be sourced)")
parser.add_argument('-s', '--set-point', type=float, required=True, help="Source value in amps or volts")
parser.add_argument('-t', '--duration', type=float, required=True, help="Number of seconds to measure for")
parser.add_argument('-c', '--compliance', type=float, required=True, help="Compliance value in amps or volts")
parser.add_argument('-n', '--nplc', type=float, default=10.0, help='NPLC value')
parser.add_argument('-f', '--front', default=False, action='store_true', help='Use the front terminals')
parser.add_argument('-w', '--four-wire', default=False, action='store_true', help='Use four wire mode')

args = parser.parse_args()

baud = 57600
terminator = bytearray.fromhex('0A').decode()

sm = k2400(visa_lib='@py', addressString=args.sm_address, terminator=terminator, serialBaud=baud, front=args.front, twoWire=not args.four_wire, quiet=True)

sm.setNPLC(args.nplc)

sm.setupDC(sourceVoltage=args.source_voltage, compliance=args.compliance, setPoint=args.set_point, senseRange='a')
sm.write(':arm:source immediate')

def measurement_callback(measurement):
    m = measurement
    # time,current,voltage,status
    print("{:0.3f},{:0.6f},{:0.6f},{:d}".format(m[2],m[1],m[0],int(m[3])),flush=True)

# print header
print('#Time [s],Current [A],Voltage [V],Status Bits',flush=True)
sm.measureUntil(t_dwell=args.duration, cb=measurement_callback)
sm.outOn(False)

