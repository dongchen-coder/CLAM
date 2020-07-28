SRC_DIR=./src
BIN_DIR=./bin

SAMPLE_TIMESTAMPED_DIR=./samples_time_stamped
SAMPLE_RANDOM_DIR=./samples_random
SAMPLE_7BIT_DIR=./lease_hardware_sampler/sampling_files/random_7bit
SAMPLE_9BIT_DIR=./lease_hardware_sampler/sampling_files/random_9bit
SAMPLE_TAG_DIR=./lease_hardware_sampler/sampling_files/random_with_tag
SAMPLE_PHASE_DIR=./lease_hardware_sampler/sampling_files/9b_phases
RESULT_DIR=./leases

TRACE_REF_SRC_DIR=./ref_trace_src
TRACE_REF_RESULT_DIR=./ref_trace_result

samples_time_stamped = 2mm_small_1000 bicg_small_1000 mvt_small_1000 3mm_small_1000 atax_small_1000 doitgen_small_1000 nussinov_small_1000

samples_random = 2mm_small_rand 3mm_small_rand atax_small_rand doitgen_small_rand mvt_small_rand nussinov_small_rand floyd_small_rand

all_bench = 2mm 3mm adi atax bicg cholesky correlation covariance deriche doitgen durbin fdtd_2d floyd_warshall gemm gemver gesummv gramschmidt heat_3d jacobi_1d jacobi_2d lu ludcmp mvt nussinov seidel_2d symm syr2d syrk trisolv trmm

samples_7bit = 2mm_small_rand 3mm_small_rand atax_small_rand doitgen_small_rand mvt_small_rand nussinov_small_rand floyd_small_rand
samples_9bit = 2mm_small_rand 3mm_small_rand atax_small_rand doitgen_small_rand mvt_small_rand nussinov_small_rand floyd_small_rand

samples_tag = 2mm_small_rand #3mm_small_rand atax_small_rand doitgen_small_rand floyd_small_rand mvt_small_rand nussinov_small_rand

phase_bench = 2mm_small_rand 3mm_small_rand atax_small_rand doitgen_small_rand floyd_small_rand mvt_small_rand nussinov_small_rand

gen:
	g++ -std=c++11 $(SRC_DIR)/rl_main.cpp -O3 -o $(BIN_DIR)/rl
	g++ -std=c++11 $(SRC_DIR)/rl_phase_main.cpp -O3 -o $(BIN_DIR)/rl_phase
	g++ -std=c++11 $(SRC_DIR)/rl_set_main.cpp -O3 -o $(BIN_DIR)/rl_set

run:
	$(foreach sample, $(samples_time_stamped), $(BIN_DIR)/rl $(SAMPLE_TIMESTAMPED_DIR)/$(sample).txt 1 > $(RESULT_DIR)/$(sample)_leases.txt ;)

run_phase:
	$(foreach sample, $(samples_time_stamped), $(BIN_DIR)/rl_phase $(SAMPLE_TIMESTAMPED_DIR)/$(sample).txt 16 > $(RESULT_DIR)/$(sample)_phase_leases.txt ;)

run_rand:
	$(foreach sample, $(samples_random), $(BIN_DIR)/rl $(SAMPLE_RANDOM_DIR)/$(sample).txt 8 > $(RESULT_DIR)/$(sample)_leases.txt ;)

run_rand_with_neg:
	$(foreach sample, $(samples_7bit), $(BIN_DIR)/rl $(SAMPLE_7BIT_DIR)/$(sample).txt 8 > $(RESULT_DIR)/$(sample)_7bit_with_neg_leases.txt ;)
	$(foreach sample, $(samples_9bit), $(BIN_DIR)/rl $(SAMPLE_9BIT_DIR)/$(sample).txt 8 > $(RESULT_DIR)/$(sample)_9bit_with_neg_leases.txt ;)

run_rand_with_set:
	$(foreach sample, $(samples_tag), $(BIN_DIR)/rl_set $(SAMPLE_TAG_DIR)/$(sample).txt 8 > $(RESULT_DIR)/$(sample)_tag_leases.txt ;)

run_rand_with_phase_set:
	$(foreach sample, $(phase_bench), $(BIN_DIR)/rl_set $(SAMPLE_PHASE_DIR)/$(sample).txt 8 > $(RESULT_DIR)/$(sample)_tag_phase_leases.txt ;)

trace_ref_run:
	$(foreach name, $(all_bench), g++ -std=c++11 $(TRACE_REF_SRC_DIR)/$(name)_ref_arr.cpp -O3 -o $(BIN_DIR)/$(name)_ref_arr_trace ;)
	$(foreach name, $(all_bench), $(BIN_DIR)/$(name)_ref_arr_trace > $(TRACE_REF_RESULT_DIR)/$(name)_ref_arr_trace_result.txt ;)
