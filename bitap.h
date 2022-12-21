 #include <string.h>
 #include <limits.h>

 const char *bitap_bitwise_search(const char *text, const char *pattern)
 {
     int m = strlen(pattern);
     unsigned long R;
     unsigned long pattern_mask[CHAR_MAX+1];
     int i;

     if (pattern[0] == '\0') return text;
     if (m > 31) return "The pattern is too long!";

     /* Initialize the bit array R */
     R = ~1;

     /* Initialize the pattern bitmasks */
     for (i=0; i <= CHAR_MAX; ++i)
         pattern_mask[i] = ~0;
     for (i=0; i < m; ++i)
         pattern_mask[pattern[i]] &= ~(1UL << i);

     for (i=0; text[i] != '\0'; ++i) {
         /* Update the bit array */
         R |= pattern_mask[text[i]];
         R <<= 1;

         if (0 == (R & (1UL << m)))
             return (text + i - m) + 1;
     }

     return NULL;
 }