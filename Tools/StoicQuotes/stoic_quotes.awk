# GRTeclyn
# Copyright 2022 The GRTL collaboration.
# Please refer to LICENSE in GRTeclyn's root directory.

function trim(value)
{
    sub(/^[[:space:]]+/, "", value)
    sub(/[[:space:]]+$/, "", value)
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

function print_box(quote, content_width, message, border, remaining,
                   split_index, box_line)
{
    content_width = 76
    message = "Stoic wisdom: " quote
    border = repeat_character("-", content_width + 2)

    print "+" border "+"
    remaining = message
    while (length(remaining) > content_width)
    {
        split_index = content_width + 1
        while (split_index > 1 && substr(remaining, split_index, 1) != " ")
        {
            --split_index
        }

        if (split_index == 1)
        {
            box_line = substr(remaining, 1, content_width)
            remaining = substr(remaining, content_width + 1)
        }
        else
        {
            box_line = substr(remaining, 1, split_index - 1)
            remaining = trim(substr(remaining, split_index + 1))
        }

        print "| " box_line \
              repeat_character(" ", content_width - length(box_line)) " |"
    }

    print "| " remaining \
          repeat_character(" ", content_width - length(remaining)) " |"
    print "+" border "+"
}

function report_error(message)
{
    print FILENAME ":" FNR ": " message > "/dev/stderr"
    invalid = 1
}

BEGIN {
    if (section != "success" && section != "failure")
    {
        print "stoic_quotes.awk: section must be success or failure" \
              > "/dev/stderr"
        invalid = 1
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

    srand()
    quote_index = 1 + int(rand() * quote_count[section])
    print_box(quotes[section, quote_index])
}
