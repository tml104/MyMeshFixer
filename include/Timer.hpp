// Copyright 2024 tml104
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <ctime>
#include <iostream>
#include <chrono>
#include <vector>

#include <spdlog/spdlog.h>

namespace Timer
{
    using MyClockType = std::chrono::time_point<std::chrono::high_resolution_clock>;

    class Timer
    {
    public:
        Timer() {Start();}

        void Start(){
            split_vec.clear();
            start_time_point = std::chrono::high_resolution_clock::now();
        }

        double GetSeconds(MyClockType c) const {
            auto duration =  std::chrono::duration_cast<std::chrono::seconds>(c - start_time_point);
            return duration.count();
        }

        void Split(std::string split_name){
            split_vec.emplace_back(split_name, std::chrono::high_resolution_clock::now());
        }

        void PrintTimes(){
            spdlog::info("===== {} =====", std::string("Timer"));

            for(auto p: split_vec){
                spdlog::info("{}: \t{:.6}", p.first, GetSeconds(p.second));
            }

            MyClockType now_time = std::chrono::high_resolution_clock::now();
            double now_sec = GetSeconds(now_time);
            spdlog::info("Total: {:.6}", now_sec);

            spdlog::info("===== {} =====", std::string("Timer End"));
        }

    private:
        MyClockType start_time_point;
        std::vector<std::pair<std::string, MyClockType>> split_vec;

    };

}