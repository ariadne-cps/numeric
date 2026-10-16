/***************************************************************************
 *            benchmark_command_line.hpp
 *
 *  Copyright  2026  Ariadne contributors
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef ARIADNE_NUMERIC_BENCHMARK_COMMAND_LINE_HPP
#define ARIADNE_NUMERIC_BENCHMARK_COMMAND_LINE_HPP

#include "utility/command_line_interface.hpp"
#include "logging/logging.hpp"

#include <algorithm>
#include <initializer_list>
#include <limits>
#include <iterator>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>

namespace Ariadne::BenchmarkCommandLine {

// Logging scheduler switches require a thread registry. The benchmarks are
// single-threaded and do not register extra application threads.
class BenchmarkThreadRegistry : public ThreadRegistryInterface {
  public:
    Bool has_threads_registered() const override { return false; }
};

inline void initialize_logging() {
    static BenchmarkThreadRegistry registry;
    auto& logger=Logger::instance();
    if (!logger.has_thread_registry_attached()) {
        logger.attach_thread_registry(&registry);
    }
    logger.configuration().set_verbosity(1u);
}


// Parse while acquiring the argument, so the standard Ariadne CLI reports
// malformed values using its normal InvalidArgumentValueException mechanism.
template<class UInt>
UInt parse_unsigned(String const& value) {
    const std::string text(value);
    if (text.empty() || text.front() == '-') {
        throw std::invalid_argument("expected nonnegative integer");
    }
    std::size_t used=0u;
    const unsigned long long parsed=std::stoull(text,&used);
    if (used!=text.size() || parsed>std::numeric_limits<UInt>::max()) {
        throw std::invalid_argument("integer out of range");
    }
    return static_cast<UInt>(parsed);
}

// The CLI itself parses only named options. Translate legacy positional
// arguments into the corresponding names without bypassing its validation,
// duplicate detection or built-in help handling.
inline List<String> normalize_legacy_arguments(
        int argc, const char* argv[],
        std::initializer_list<std::string_view> positional_options,
        std::initializer_list<std::string_view> valued_options)
{
    List<String> arguments;
    arguments.append(String(argv[0]));
    SizeType positional_index=0u;

    for (int i=1; i<argc; ++i) {
        const String argument(argv[i]);
        const std::string_view view(argument);
        if (std::find(valued_options.begin(),valued_options.end(),view)
                != valued_options.end()) {
            arguments.append(argument);
            if (i+1<argc) {
                arguments.append(String(argv[++i]));
            }
        } else if (!argument.empty() && argument.front()!='-'
                   && positional_index<positional_options.size()) {
            auto position=positional_options.begin();
            std::advance(position,static_cast<std::ptrdiff_t>(positional_index));
            arguments.append(String("--")+String(std::string(*position)));
            arguments.append(argument);
            ++positional_index;
        } else {
            arguments.append(argument);
        }
    }
    return arguments;
}

inline bool requests_help(List<String> const& arguments) {
    return std::find(arguments.begin(),arguments.end(),String("--help"))!=arguments.end()
        || std::find(arguments.begin(),arguments.end(),String("-h"))!=arguments.end();
}

} // namespace Ariadne::BenchmarkCommandLine

#endif // ARIADNE_NUMERIC_BENCHMARK_COMMAND_LINE_HPP
