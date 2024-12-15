//
//  logger_V2.hpp
//  Final_Project_V3
//
//  Created by jonathon hope on 12/14/24.
//

#pragma once

#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <mutex>
#include <string>

using namespace std;

// Enum to represent log levels
enum LogLevel
{
    HIT_POINCARE,
    HIT_MOON,
    POS_VS_NEG,
    MEASUREMENTS,
    SC_STATE,
    PO_DATA,
    MAX_DELTA_V,
    TRANSFER_INFO,
    BRANCH_COUNT,
    FUNCTION_START,
    FUNCTION_END,
    PROGRAM_START,
    PROGRAM_END
};

class Logger
{
private:
    // Static pointer to the Singleton instance
    static Logger *instancePtr;

    // Mutex to ensure thread safety
    static mutex mtx;

    ofstream logFile; // File stream for the log file

    // Private constructor
    Logger(const string &filename);

    // Private method to open the log file
    void openLogFile(const string &filename);

public:
    // Deleting the copy constructor to prevent copies
    Logger(const Logger &obj) = delete;

    // Static method to get the Singleton instance
    static Logger *getInstance(const string &filename = "");

    // Destructor: Closes the log file
    ~Logger();

    // Set log file for the logger
    void setLogFile(const string &filename);

    // Logs a message with a given log level
    void log(LogLevel level, const string &message);

    // Converts log level to a string for output
    string levelToString(LogLevel level);
};
