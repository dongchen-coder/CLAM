SRC_DIR=./src
BIN_DIR=./bin
SAMPLE_TIMESTAMPED_DIR=./samples_time_stamped
RESULT_DIR=./leases

samples_time_stamped = 2mm_small_1000 bicg_small_1000 mvt_small_1000 3mm_small_1000 atax_small_1000 doitgen_small_1000 nussinov_small_1000

gen:
	g++ -std=c++11 $(SRC_DIR)/rl_main.cpp -O3 -o $(BIN_DIR)/rl
	g++ -std=c++11 $(SRC_DIR)/rl_phase_main.cpp -O3 -o $(BIN_DIR)/rl_phase

run:
	$(foreach sample, $(samples_time_stamped), $(BIN_DIR)/rl $(SAMPLE_TIMESTAMPED_DIR)/$(sample).txt 1 > $(RESULT_DIR)/$(sample)_leases.txt ;)

run_phase:
	$(foreach sample, $(samples_time_stamped), $(BIN_DIR)/rl_phase $(SAMPLE_TIMESTAMPED_DIR)/$(sample).txt 16 > $(RESULT_DIR)/$(sample)_phase_leases.txt ;)
