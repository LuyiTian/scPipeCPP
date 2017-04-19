#include "timing.h"

using namespace std::chrono;

void Timer::start() {
	start_time = steady_clock::now();
}

std::string Timer::end() {
	end_time = steady_clock::now();
	std::string elapsed_string;
	
	if (time_mode == "seconds") {
		seconds elapsed_time = duration_cast<seconds>(end_time - start_time);
		elapsed_string = std::to_string(elapsed_time.count()) + " seconds";
	} else if (time_mode == "minutes") {
		minutes elapsed_time = duration_cast<minutes>(end_time - start_time);
		elapsed_string = std::to_string(elapsed_time.count()) + " minutes";
	} else if (time_mode == "hours") {
		hours elapsed_time = duration_cast<hours>(end_time - start_time);
		elapsed_string = std::to_string(elapsed_time.count()) + " hours";
	}

	return elapsed_string;
}

void Timer::mode(std::string mode) {
	time_mode = mode;
}
