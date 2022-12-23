 #include <string.h>
 #include <limits.h>
 #include <stdint.h>

 const char *bitap_bitwise_search(const char *haystack, size_t haystack_length,
                                  const char *needle, size_t needle_length)
{
    size_t R;
    size_t pattern_mask[CHAR_MAX+1];
    size_t i;

    if (needle_length == 0) return haystack;
    if (needle_length > (sizeof(size_t) * 8 -1 ))
        return "The pattern is too long!";

    /* Initialize the bit array R */
    R = ~1;

    /* Initialize the pattern bitmasks */
    memset(pattern_mask, 0xff, sizeof(size_t) * CHAR_MAX);
    for (i=0; i < needle_length; ++i)
        pattern_mask[needle[(uint8_t)i]] &= ~(1UL << i);

    for (i=0; i < haystack_length; ++i) {
        /* Update the bit array */
        R |= pattern_mask[haystack[(uint8_t)i]];
        R <<= 1;

        if (0 == (R & (1UL << needle_length)))
            return (haystack + i - needle_length) + 1;
     }

    return NULL;
 }