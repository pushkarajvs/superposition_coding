all: single_user_12rates_far_op single_user_12rates_near_op superposition_coding_op

single_user_12rates_far_op: single_user_12rates_with_channel_estimate.cpp
	g++ -o single_user_12rates_far_op -DFAR=1 single_user_12rates_with_channel_estimate.cpp -litpp

single_user_12rates_near_op: single_user_12rates_with_channel_estimate.cpp
	g++ -o single_user_12rates_near_op -DNEAR=1 single_user_12rates_with_channel_estimate.cpp -litpp

superposition_coding_op: superposition_coding_for_multipath_channel_estimate.cpp
	g++ -o superposition_coding_op superposition_coding_for_multipath_channel_estimate.cpp -litpp

.PHONY: clean
clean:
	$(RM) superposition_coding_op single_user_12rates_far_op single_user_12rates_near_op
