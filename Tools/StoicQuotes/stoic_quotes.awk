# GRTeclyn
# Copyright 2022 The GRTL collaboration.
# Please refer to LICENSE in GRTeclyn's root directory.

function trim(value)
{
    sub(/^[[:space:]]+/, "", value)
    sub(/[[:space:]]+$/, "", value)
    return value
}

function cpp_escape(value)
{
    gsub(/\\/, "\\\\", value)
    gsub(/"/, "\\\"", value)
    return value
}

function repeat_character(character, count, result, character_index)
{
    result = ""
    for (character_index = 0; character_index < count; ++character_index)
    {
        result = result character
    }
    return result
}

function report_error(message)
{
    print FILENAME ":" FNR ": " message > "/dev/stderr"
    invalid = 1
}

BEGIN {
    if (mode != "generate" && mode != "random")
    {
        print "stoic_quotes.awk: mode must be generate or random" > "/dev/stderr"
        exit 1
    }
}

{
    line = trim($0)

    if (line == "" || line ~ /^#/)
    {
        next
    }

    if (line ~ /^\[[^]]+\]$/)
    {
        current_section = tolower(substr(line, 2, length(line) - 2))
        if (current_section != "success" && current_section != "failure")
        {
            report_error("unknown section " line)
            current_section = ""
        }
        next
    }

    if (current_section == "")
    {
        report_error("quote appears before a valid section")
        next
    }

    quotes[current_section, ++quote_count[current_section]] = line
}

END {
    if (quote_count["success"] == 0)
    {
        report_error("empty [success] section")
    }
    if (quote_count["failure"] == 0)
    {
        report_error("empty [failure] section")
    }
    if (invalid)
    {
        exit 1
    }

    if (mode == "random")
    {
        srand()
        quote_index = 1 + int(rand() * quote_count[section])
        message = " Stoic wisdom: " quotes[section, quote_index] " "
        border = repeat_character("-", length(message))
        print "+" border "+"
        print "|" message "|"
        print "+" border "+"
        exit
    }

    print "// Generated from StoicQuotes.txt. Do not edit directly."
    print "#include \"StoicQuotes.hpp\""
    print ""
    print "#include <array>"
    print "#include <cstddef>"
    print "#include <random>"
    print "#include <string_view>"
    print ""
    print "namespace StoicQuotes"
    print "{"
    print "namespace"
    print "{"

    section_names[1] = "success"
    section_names[2] = "failure"
    for (section_index = 1; section_index <= 2; ++section_index)
    {
        current_section = section_names[section_index]
        print "inline constexpr std::array " current_section " = {"
        for (quote_index = 1;
             quote_index <= quote_count[current_section]; ++quote_index)
        {
            suffix = quote_index < quote_count[current_section] ? "," : "};"
            print "    std::string_view{\"" \
                  cpp_escape(quotes[current_section, quote_index]) "\"}" suffix
        }
        print ""
    }

    print "} // namespace"
    print ""
    print "std::string_view random_quote(bool succeeded)"
    print "{"
    print "    static std::mt19937 generator("
    print "        static_cast<std::mt19937::result_type>(std::random_device{}()));"
    print ""
    print "    if (succeeded)"
    print "    {"
    print "        std::uniform_int_distribution<std::size_t> distribution("
    print "            0, success.size() - 1);"
    print "        return success[distribution(generator)];"
    print "    }"
    print ""
    print "    std::uniform_int_distribution<std::size_t> distribution("
    print "        0, failure.size() - 1);"
    print "    return failure[distribution(generator)];"
    print "}"
    print ""
    print "} // namespace StoicQuotes"
}
