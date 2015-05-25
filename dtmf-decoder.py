'''
A python implementation of the Goertzel algorithm to decode DTMF tones.
The wave file is split into bins and each bin is analyzed
for all the DTMF frequencies. The method run() will return a numeric
representation of the DTMF tone.
'''
import wave
import struct
import sys
import math

ITU_Q23_ROWS = ( 697.0,  770.0,  852.0,  941.0) # Hz
ITU_Q23_COLS = (1209.0, 1336.0, 1477.0, 1633.0) # Hz
ITU_Q23_FREQS = ITU_Q23_ROWS + ITU_Q23_COLS

ITU_Q23_KEYS = {# Row    Col     Key
                (697.0, 1209.0): '1',
                (697.0, 1336.0): '2',
                (697.0, 1477.0): '3',
                (697.0, 1633.0): 'A',
                (770.0, 1209.0): '4',
                (770.0, 1336.0): '5',
                (770.0, 1477.0): '6',
                (770.0, 1633.0): 'B',
                (852.0, 1209.0): '7',
                (852.0, 1336.0): '8',
                (852.0, 1477.0): '9',
                (852.0, 1633.0): 'C',
                (941.0, 1209.0): '*',
                (941.0, 1336.0): '0',
                (941.0, 1477.0): '#',
                (941.0, 1633.0): 'D',
                }


class pygoertzel_dtmf:
    def __init__(self, samplerate, freq=ITU_Q23_FREQS):
        self.samplerate = samplerate
        self.goertzel_freq = freq
        self.s_prev = {}
        self.s_prev2 = {}
        self.totalpower = {}
        self.N = {}
        self.coeff = {}
        # create goertzel parameters for each frequency so that
        # all the frequencies are analyzed in parallel
        for k in self.goertzel_freq:
            self.s_prev[k] = 0.0
            self.s_prev2[k] = 0.0
            self.totalpower[k] = 0.0
            self.N[k] = 0.0
            normalizedfreq = k / self.samplerate
            self.coeff[k] = 2.0*math.cos(2.0 * math.pi * normalizedfreq)

    def strongest_freq(self, freqs, candidate_freqs):
        freq = 0.0
        freq_v = 0.0
        for f in candidate_freqs:
            if freqs[f] > freq_v:
                freq_v = freqs[f]
                freq = f
        return freq

    def __get_number(self, freqs):
        hifreq = self.strongest_freq(freqs, reversed(ITU_Q23_COLS))
        lofreq = self.strongest_freq(freqs, reversed(ITU_Q23_ROWS))

        if lofreq and hifreq:
            return ITU_Q23_KEYS[(lofreq, hifreq)]

    def run(self, sample):
        freqs = {}
        for freq in self.goertzel_freq:
            s = sample + (self.coeff[freq] * self.s_prev[freq]) - self.s_prev2[freq]
            self.s_prev2[freq] = self.s_prev[freq]
            self.s_prev[freq] = s
            self.N[freq]+=1
            power = (self.s_prev2[freq]*self.s_prev2[freq]) + (self.s_prev[freq]*self.s_prev[freq]) - (self.coeff[freq]*self.s_prev[freq]*self.s_prev2[freq])
            self.totalpower[freq]+=sample*sample
            if (self.totalpower[freq] == 0):
                self.totalpower[freq] = 1
            freqs[freq] = power / self.totalpower[freq] / self.N[freq]
        return self.__get_number(freqs)

if __name__ == '__main__':
    # load wav file
    wav = wave.open(sys.argv[1], 'r')
    (nchannels, sampwidth, framerate, nframes, comptype, compname) = wav.getparams()
    frames = wav.readframes(nframes * nchannels)
    # convert wave file to array of integers
    frames = struct.unpack_from("%dH" % nframes * nchannels, frames)
    # if stereo get left/right
    if nchannels == 2:
        left = [frames[i] for i in range(0,len(frames),2)]
        right = [frames[i] for i in range(1,len(frames),2)]
    else:
        left = frames
        right = left
    binsize = 400
    # Split the bin in 4 to average out errors due to noise
    binsize_split = 4
    prevvalue = ""
    prevcounter = 0
    for i in range(0,len(left)-binsize,binsize/binsize_split):
        goertzel = pygoertzel_dtmf(framerate)
        for j in left[i:i+binsize]:
            value = goertzel.run(j)
        if value==prevvalue:
            prevcounter+=1
            if prevcounter==10:
                print value
        else:
            prevcounter=0
            prevvalue=value
