# -*- coding: utf-8 -*-
from bitarray import bitarray
from bitarray_utils import intToBin, binToInt, grayToBin, grayToInt

def nrange(start, end):
	r = [x for x in range(start, end) ]
	r.reverse()
	return r


GlobalConfigAfterReset = bitarray('1100011110010111111101011010110111001011011110001101000100110011101110110110011001011111001110111110111011111000111111111110010111101111111001111111011000010011111100001101000000010110')

## Contains parameters and methods related to the global operation of one ASIC. 
class AsicGlobalConfig(bitarray):
	## Constructor. 
	# Defines and sets all fields to default values. Most important fields are:
	# 
	def __init__(self, initial=None, endian="big"):
		super(AsicGlobalConfig, self).__init__()

		self.__fields = {
			"r_clk_en"	 : [2, 1, 0],
			"tac_refresh_en"	 : [3],
			"tac_refresh_period"	 : [6, 5, 7],
			"counter_en"	 : [8],
			"counter_period"	 : [12, 11, 10, 9],
			"tx_is_debug"	 : [13],
			"tx_ddr"	 : [14],
			"tx_nlinks"	 : [16, 15],
			"tx_mode"	 : [18, 17],
			"veto_en"	 : [19],
			"fe_test_pulse"	 : [20],
			"tgr_enable"	 : [21],
			"tgr_wnd_a"	 : [33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22],
			"tgr_wnd_z"	 : [45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34],
			"tgr_opportunistic_tx"	 : [46],
			"tp_ext_en"	 : [47],
			"tp_int_en"	 : [48],
			"tp_int_inv"	 : [49],
			"tp_int_period"	 : [59, 58, 57, 56, 55, 54, 53, 52, 51, 50],
			"tp_int_length"	 : [69, 68, 67, 66, 65, 64, 63, 62, 61, 60],
			"vb1_bias_ig"	 : [75, 74, 73, 72, 71, 70],
			"main_global_dac"	 : [80, 79, 78, 77, 76],
			"vdiscv_bias_ig"	 : [84, 83, 82, 81],
			"tdc_global_dac"	 : [90, 89, 88, 87, 86, 85],
			"sar_vbias"	 : [96, 95, 94, 93, 92, 91],
			"iref_cs"	 : [101, 100, 99, 98, 97],
			"tdc_vcasn"	 : [105, 104, 103, 102],
			"fe_ib2"	 : [111, 110, 109, 108, 107, 106],
			"disc_sf_bias"	 : [117, 116, 115, 114, 113, 112],
			"vdisci_bias_ig"	 : [123, 122, 121, 120, 119, 118],
			"debug_buff_bias"	 : [126, 125, 124],
			"vdifffoldcas_dac_ig"	 : [132, 131, 130, 129, 128, 127],
			"imirror_bias_top"	 : [137, 136, 135, 134, 133],
			"ilsb"	 : [138, 139, 140, 141, 142],
			"tac_vcasp"	 : [143, 144, 145, 146, 147],
			"v_cal_tp_top"	 : [152, 151, 150, 149, 148],
			"v_cal_diff_bias_ig"	 : [157, 156, 155, 154, 153],
			"disc_lsb_t1"	 : [163, 162, 161, 160, 159, 158],
			"disc_lsb_t2"	 : [169, 168, 167, 166, 165, 164],
			"v_att_diff_bias_ig"	 : [175, 174, 173, 172, 171, 170],
			"v_ref_diff_bias_ig"	 : [181, 180, 179, 178, 177, 176],
			"v_cal_ref_ig"	 : [186, 185, 184, 183, 182],
			"v_ipostamp_top"	 : [191, 190, 189, 188, 187],
			"v_ipostamp_e_top"	 : [196, 195, 194, 193, 192],
			"disc_lsb_e"	 : [202, 201, 200, 199, 198, 197],
			"v_att_ls_bias_ig"	 : [208, 207, 206, 205, 204, 203],
			"v_integ_ref_ig"	 : [214, 213, 212, 211, 210, 209],
			"cg2_sw_n"	 : [215],
			"fetp_en"	 : [220],
			"input_pol"	 : [221]
		}

		if initial is not None:
			# We have an initial value, let's just go with it!
			self[0:222] = bitarray(initial)	;

			return None


		self[0:222] = bitarray('100000111101110111111100110100101111111111100110001110110111101011111111111011011001101111101101110100111010010000000111100101100001101001110101111101110000000000000000000000000000000000000000000000000001000000000000000100')
		
		# WARNING These values were set by experimental adjustment
		# tdc global dac adjust due to mismatch
		self.setValue("tdc_global_dac", 42)
		# main global dac adjustment due to mismatch
		self.setValue("main_global_dac", 0b10111) ## WARNING

		# INFO: values set from simulation
		self.setValue("main_global_dac", 0b10111)
		self.setValue("tdc_global_dac", 0b110100)
		
		# WARNING: Pushing these values to maximum seems to lead to better SPTR
		# but also push power consumption. Needs more study
		#self.setValue("fe_ib2", 0b0100000)
		#self.setValue("fe_ib2", 0b0000000) ## WARNING
		#self.setValue("disc_sf_bias", 0)
		
		
		# WARNING These seem to be a reasonable compromise on having widest range for these discriminators
		# but was obtained with a small sampe
		self.setValue("disc_lsb_t2", 48)
		self.setValue("disc_lsb_e", 40)

		# Default FETP settings
		self.setValue("v_cal_tp_top", 1)
		self.setValue("v_cal_diff_bias_ig", 0)
		self.setValue("v_cal_ref_ig", 31);

		# Disable the counter and set it to a reasonably long period
		self.setValue("counter_en", 0b0)
		self.setValue("counter_period", 0b110)

		return None

	def __deepcopy__(self, memo):
		return AsicGlobalConfig(initial=self)
		

	## Set the value of a given parameter as an integer
	# @param key  String with the name of the parameter to be set
	# @param value  Integer corresponding to the value to be set	
	def setValue(self, key, value):
		b = intToBin(value, len(self.__fields[key]))
		self.setBits(key, b)

	## Set the value of a given parameter as a bitarray
	# @param key  String with the name of the parameter to be set
	# @param value  Bitarray corresponding to the value to be set		
	def setBits(self, key, value):
		index = self.__fields[key]
		assert len(value) == len(index)
		for a,b in enumerate(index):
			self[221 - b] = value[a]

	## Returns the value of a given parameter as a bitarray
	# @param key  String with the name of the parameter to be returned	
	def getBits(self, key):
		index = self.__fields[key]
		value = bitarray(len(index))
		for a,b in enumerate(index):
			value[a] = self[221 - b]
		return value

	## Returns the value of a given parameter as an integer
	# @param key  String with the name of the parameter to be returned	
	def getValue(self, key):
		return binToInt(self.getBits(key))

	## Prints the content of all parameters as a bitarray	
	def printAllBits(self):
		for key in self.__fields.keys():
			print key, " : ", self.getBits(key)

	## Prints the content of all parameters as integers
	def printAllValues(self):
		unsorted = [ (min(bitList), name) for name, bitList in self.__fields.items() ]
		unsorted.sort()
		for b, key in unsorted:
			bitList = self.__fields[key]
			l = bitList[0]
			r = bitList[-1]
			print "%30s : %3d : %20s : %d..%d" % (key, self.getValue(key), self.getBits(key), l, r)

	## Returns all the keys (variables) in this class
	def getKeys(self):
		return self.__fields.keys()

## Contains parameters and methods related to the operation of one channel of the ASIC. 
class AsicChannelConfig(bitarray):
	## Constructor
	# Defines and sets all fields to default values. Most important fields are:
	# 
	def __init__(self, initial=None, endian="big"):
		super(AsicChannelConfig, self).__init__()

		self.__fields = {
			"metastability"	 : [1, 0],
			"trigger_mode_1"	 : [3, 2],
			"trigger_mode_2_t1"	 : [5, 4],
			"trigger_mode_2_q"	 : [7, 6],
			"trigger_mode_2_t2"	 : [9, 8],
			"trigger_mode_2_v"	 : [11, 10],
			"trigger_mode_2_b"	 : [14, 13, 12],
			"debug"	 : [16, 15],
			"branch_en"	 : [19, 18, 17],
			"min_intg_time"	 : [26, 25, 24, 23, 22, 21, 20],
			"max_intg_time"	 : [33, 32, 31, 30, 29, 28, 27],
			"trigger_b_latched"	 : [34],
			"counter_mode"	 : [38, 37, 36, 35],
			"dead_time"	 : [44, 43, 42, 41, 40, 39],
			"min_valid_time"	 : [48, 47, 46, 45],
			"max_valid_time"	 : [52, 51, 50, 49],
			"max_valid_time_en"	 : [53],
			"tac_max_age"	 : [58, 57, 56, 55, 54],
			"reset_n2"	 : [59],
			"stable2"	 : [60],
			"qtx2"	 : [61],
			"conv2"	 : [62],
			"tac_rd_ext"	 : [63],
			"fe_delay"	 : [68, 116, 114, 64, 66 ],
			"vth_e"	 : [73, 71, 69, 67, 65, 83],
			"hytst_t1"	 : [74, 72, 70],
			"vth_t2"	 : [81, 79, 77, 75, 93, 91],
			"hytst_t2"	 : [80, 78, 76],
			"hyst_e"	 : [86, 84, 82],
			"vth_t1"	 : [89, 87, 85, 103, 101, 99],
			"sh_e"	 : [90, 88],
			"gain_e"	 : [92, 94],
			"baseline_t"	 : [97, 95, 111, 110, 108, 106],
			"baseline_e"	 : [96, 98, 100],
			"gain_t"	 : [102, 105],
			"output_eh_h"	 : [104],
			"integ_source_sw"	 : [107, 109],
			"fetp_en"	 : [126, 112],
			"att"	 : [113, 115, 117],
			"hyst_en_n"	 : [118],
			"qdc_en"	 : [119],
			"tdc_t_iref"	 : [120, 121, 122, 123],
			"tdc_e_iref"	 : [124, 125, 128, 129],
			"input_pol"	 : [127]
		}

		if initial is not None:
			# We have an initial value, let's just go with it!
			self[0:130] = bitarray(initial)
			return None

		# Specify default value
		self[0:130] = bitarray('0010000000110010101110010000010111101010101011101001101001100100100000001010100000000000000010000000100000010011100101100101110000')

		# Disable shaping by default
		# Disable shaping by default
		self.setValue("sh_e", 0b00)
		
		# Default integration windows: fixed ~300 ns
		self.setValue("min_intg_time", 34)
		self.setValue("max_intg_time", 34)

		# Avoid powers of 2
		self.setValue("tac_max_age", 30)
		
		# This setting gives better time resolution
		self.setValue("fe_delay", 14)

		# This setting has better linearity
		self.setValue("att", 1)
		
		# Disable max valid time
		# It should only be used with very short signals
		self.setValue("max_valid_time_en", 0b0)
		# Trigger_T2 on falling edge of dead time does not work well
		self.setValue("dead_time", 1)
		
		# Enable 3 SAR conversions
		self.setValue("branch_en", 0b111)

		return None

	def __deepcopy__(self, memo):
		return AsicChannelConfig(initial=self)

	## Set the value of a given parameter as an integer
	# @param key  String with the name of the parameter to be set
	# @param value  Integer corresponding to the value to be set		
	def setValue(self, key, value):
		b = intToBin(value, len(self.__fields[key]))
		self.setBits(key, b)

	## Set the value of a given parameter as a bitarray
	# @param key  String with the name of the parameter to be set
	# @param value  Bitarray corresponding to the value to be set		
	def setBits(self, key, value):
		index = self.__fields[key]
		assert len(value) == len(index)
		for a,b in enumerate(index):
			self[129 - b] = value[a]

	## Returns the value of a given parameter as a bitarray
	# @param key  String with the name of the parameter to be returned	
	def getBits(self, key):
		index = self.__fields[key]
		value = bitarray(len(index))
		for a,b in enumerate(index):
			value[a] = self[129 - b]
		return value

	## Returns the value of a given parameter as an integer
	# @param key  String with the name of the parameter to be returned	
	def getValue(self, key):
		return binToInt(self.getBits(key))

	# Prints the content of all parameters as a bitarray		
	def printAllBits(self):
		for key in self.__fields.keys():
			print key, " : ", self.getBits(key)
		
	## Prints the content of all parameters as integers
	def printAllValues(self):
		unsorted = [ (min(bitList), name) for name, bitList in self.__fields.items() ]
		unsorted.sort()
		for b, key in unsorted:
			bitList = self.__fields[key]
			l = bitList[0]
			r = bitList[-1]
			print "%30s : %3d : %20s : %d..%d" % (key, self.getValue(key), self.getBits(key), l, r)

	## Set the baseline value in units of ADC (63 to 0)
	def setBaseline(self, v):
		self.__baseline = v

	## Returns the baseline value for this channel
	def getBaseline(self):
		return self.__baseline

	## Returns all the keys (variables) in this class
	def getKeys(self):
		return self.__fields.keys()

## A class containing instances of AsicGlobalConfig and AsicChannelConfig
#, as well as 2 other bitarrays related to test pulse configuration. Is related to one given ASIC.
class AsicConfig:
	def __init__(self):
		self.channelConfig = [ AsicChannelConfig() for x in range(16) ]
		self.globalConfig = AsicGlobalConfig()
		return None



class ConfigurationError:
	pass

class ConfigurationErrorBadAck(ConfigurationError):
	def __init__(self, portID, slaveID, asicID, value):
		self.addr = (value, portID, slaveID, asicID)
		self.errType = value
	def __str__(self):
		return "Bad ACK (%d) when configuring ASIC at port %2d, slave %2d, asic %2d"  % self.addr

class ConfigurationErrorBadCRC(ConfigurationError):
	def __init__(self, portID, slaveID, asicID):
		self.addr = (portID, slaveID, asicID)
	def __str__(self):
		return "Received configuration datta with bad CRC from ASIC at port %2d, slave %2d, asic %2d" % self.addr

class ConfigurationErrorStuckHigh(ConfigurationError):
	def __init__(self, portID, slaveID, asicID):
		self.addr = (portID, slaveID, asicID)
	def __str__(self):
		return "MOSI stuck high from ASIC at port %2d, slave %2d, asic %2d" % self.addr

class ConfigurationErrorGeneric(ConfigurationError):
	def __init__(self, portID, slaveID, asicID, value):
		self.addr = (value, portID, slaveID, asicID)
	def __str__(self):
		return "Unexpected configuration error %02X from ASIC at port %2d, slave %2d, asic %2d" % self.addr

class ConfigurationErrorBadRead(ConfigurationError):
	def __init__(self, portID, slaveID, asicID, written, read):
		self.data = (portID, slaveID, asicID, written, read)
	def __str__(self):
		return "Configuration readback failed for ASIC at port %2d, slave %2d, asic %2d: wrote %s, read %s" % self.data

class ConfigurationErrorBadReply(ConfigurationError):
	def __init__(self, expected, actual):
		self.data = (expected, actual)
	def __str__(self):
		return "Bad reply for ASIC configuration command: expected %d bytes, got %d bytes" % self.data


