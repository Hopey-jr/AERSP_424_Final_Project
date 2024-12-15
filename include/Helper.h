//
//  Helper.h
//  Final_Project_V3
//
//  Created by jonathon hope on 12/14/24.
//

#pragma once
#define __FILENAME__ (__builtin_strrchr(__FILE__, '/') ? __builtin_strrchr(__FILE__, '/') + 1 : __FILE__)

#define LOG(fmt, ...) printf("%s:%d | " fmt, __FILENAME__, __LINE__, __VA_ARGS__);
