# Summarize Context

Summarize the current conversation context.

## Usage
```
/summarize <num-words> <instructions>
```

**Parameters:**
- `num-words` - Target word count for the summary
- `instructions` - Focus area or specific aspects to emphasize in the summary

## Examples

```bash
# Summarize in 100 words focusing on technical decisions
/summarize 100 technical decisions made

# Summarize in 50 words focusing on errors encountered
/summarize 50 errors and fixes

# Summarize in 200 words focusing on the current task
/summarize 200 what we're currently working on
```

## Instructions

When this command is invoked with arguments `$ARGUMENTS`:

1. Parse the arguments:
   - First argument is the target word count (number)
   - Remaining arguments are the focus instructions

2. Review the entire conversation context including:
   - User requests and goals
   - Actions taken
   - Files modified
   - Errors encountered
   - Current state and pending work

3. Generate a summary that:
   - Is approximately the specified word count (within 10%)
   - Focuses on the aspects specified in the instructions
   - Is written in clear, concise prose
   - Highlights the most important information relevant to the instructions

4. Output the summary directly to the user.
