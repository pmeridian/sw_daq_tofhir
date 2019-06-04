TARGET=$1

mkdir ${TARGET}
#cp data/config_common.ini ${TARGET}/config.ini
./make_bias_calibration_table -o ${TARGET}/bias_calibration.tsv
./make_simple_bias_settings_table \
	--config ${TARGET}/config.ini \
	--prebd 40 \
	--bd 65.33 \
	--over 4.0 \
	--offset 0.0 \
	-o ${TARGET}/bias_settings.tsv 

./make_fake_disc_calibration_table -o ${TARGET}/disc_calibration.tsv
./make_simple_disc_settings_table \
	--config ${TARGET}/config.ini \
	--vth_t1 20 \
	--vth_t2 20 \
	--vth_e 20 \
	-o ${TARGET}/disc_settings.tsv

./make_simple_channel_map \
	--port 0 \
	--slave 0 \
	--trigger-port 0 \
	-o ${TARGET}/map

./acquire_tdc_calibration \
	--config ${TARGET}/config.ini \
	-o ${TARGET}/tdc_calibration

./process_tdc_calibration \
	--config ${TARGET}/config.ini \
	-i ${TARGET}/tdc_calibration \
	-o ${TARGET}/tdc_calibration



