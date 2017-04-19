#pragma once
#include <string>
#include <chrono>

class Timer {
public:
	void start();
	std::string end();
	void mode(std::string mode); // mode one of "seconds", "minutes", "hours"
private:
	std::chrono::time_point<std::chrono::steady_clock> start_time, end_time;
	std::string time_mode = "seconds";
};