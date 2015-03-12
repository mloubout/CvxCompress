#include <string.h>

/*
 * This code converts an ASCII string into an EBCDIC Code Page 500 string.
 */

void Convert_ASCII_To_EBCDIC(
	char* ASCII_str,
	char* EBCDIC_str,
	int num_chars
	)
{
	int i;
	for (i = 0;  i < num_chars;  ++i)
	{
		switch (ASCII_str[i])
		{
		case '[':
			EBCDIC_str[i] = 74;
			break;
		case '.':
			EBCDIC_str[i] = 75;
			break;
		case '<':
			EBCDIC_str[i] = 76;
			break;
		case '(':
			EBCDIC_str[i] = 77;
			break;
		case '+':
			EBCDIC_str[i] = 78;
			break;
		case '!':
			EBCDIC_str[i] = 79;
			break;
		case '&':
			EBCDIC_str[i] = 80;
			break;
		case ']':
			EBCDIC_str[i] = 90;
			break;
		case '$':
			EBCDIC_str[i] = 91;
			break;
		case '*':
			EBCDIC_str[i] = 92;
			break;
		case ')':
			EBCDIC_str[i] = 93;
			break;
		case ';':
			EBCDIC_str[i] = 94;
			break;
		case '^':
			EBCDIC_str[i] = 95;
			break;
		case '-':
			EBCDIC_str[i] = 96;
			break;
		case '/':
			EBCDIC_str[i] = 97;
			break;
		case '|':
			EBCDIC_str[i] = 106;
			break;
		case ',':
			EBCDIC_str[i] = 107;
			break;
		case '%':
			EBCDIC_str[i] = 108;
			break;
		case '_':
			EBCDIC_str[i] = 109;
			break;
		case '>':
			EBCDIC_str[i] = 110;
			break;
		case '?':
			EBCDIC_str[i] = 111;
			break;
		case '`':
			EBCDIC_str[i] = 121;
			break;
		case ':':
			EBCDIC_str[i] = 122;
			break;
		case '#':
			EBCDIC_str[i] = 123;
			break;
		case '@':
			EBCDIC_str[i] = 124;
			break;
		case '\'':
			EBCDIC_str[i] = 125;
			break;
		case '=':
			EBCDIC_str[i] = 126;
			break;
		case '"':
			EBCDIC_str[i] = 127;
			break;
		case 'a':
			EBCDIC_str[i] = 129;
			break;
		case 'b':
			EBCDIC_str[i] = 130;
			break;
		case 'c':
			EBCDIC_str[i] = 131;
			break;
		case 'd':
			EBCDIC_str[i] = 132;
			break;
		case 'e':
			EBCDIC_str[i] = 133;
			break;
		case 'f':
			EBCDIC_str[i] = 134;
			break;
		case 'g':
			EBCDIC_str[i] = 135;
			break;
		case 'h':
			EBCDIC_str[i] = 136;
			break;
		case 'i':
			EBCDIC_str[i] = 137;
			break;
		case 'j':
			EBCDIC_str[i] = 145;
			break;
		case 'k':
			EBCDIC_str[i] = 146;
			break;
		case 'l':
			EBCDIC_str[i] = 147;
			break;
		case 'm':
			EBCDIC_str[i] = 148;
			break;
		case 'n':
			EBCDIC_str[i] = 149;
			break;
		case 'o':
			EBCDIC_str[i] = 150;
			break;
		case 'p':
			EBCDIC_str[i] = 151;
			break;
		case 'q':
			EBCDIC_str[i] = 152;
			break;
		case 'r':
			EBCDIC_str[i] = 153;
			break;
		case '~':
			EBCDIC_str[i] = 161;
			break;
		case 's':
			EBCDIC_str[i] = 162;
			break;
		case 't':
			EBCDIC_str[i] = 163;
			break;
		case 'u':
			EBCDIC_str[i] = 164;
			break;
		case 'v':
			EBCDIC_str[i] = 165;
			break;
		case 'w':
			EBCDIC_str[i] = 166;
			break;
		case 'x':
			EBCDIC_str[i] = 167;
			break;
		case 'y':
			EBCDIC_str[i] = 168;
			break;
		case 'z':
			EBCDIC_str[i] = 169;
			break;
		case '{':
			EBCDIC_str[i] = 192;
			break;
		case 'A':
			EBCDIC_str[i] = 193;
			break;
		case 'B':
			EBCDIC_str[i] = 194;
			break;
		case 'C':
			EBCDIC_str[i] = 195;
			break;
		case 'D':
			EBCDIC_str[i] = 196;
			break;
		case 'E':
			EBCDIC_str[i] = 197;
			break;
		case 'F':
			EBCDIC_str[i] = 198;
			break;
		case 'G':
			EBCDIC_str[i] = 199;
			break;
		case 'H':
			EBCDIC_str[i] = 200;
			break;
		case 'I':
			EBCDIC_str[i] = 201;
			break;
		case '}':
			EBCDIC_str[i] = 208;
			break;
		case 'J':
			EBCDIC_str[i] = 209;
			break;
		case 'K':
			EBCDIC_str[i] = 210;
			break;
		case 'L':
			EBCDIC_str[i] = 211;
			break;
		case 'M':
			EBCDIC_str[i] = 212;
			break;
		case 'N':
			EBCDIC_str[i] = 213;
			break;
		case 'O':
			EBCDIC_str[i] = 214;
			break;
		case 'P':
			EBCDIC_str[i] = 215;
			break;
		case 'Q':
			EBCDIC_str[i] = 216;
			break;
		case 'R':
			EBCDIC_str[i] = 217;
			break;
		case '\\':
			EBCDIC_str[i] = 224;
			break;
		case 'S':
			EBCDIC_str[i] = 226;
			break;
		case 'T':
			EBCDIC_str[i] = 227;
			break;
		case 'U':
			EBCDIC_str[i] = 228;
			break;
		case 'V':
			EBCDIC_str[i] = 229;
			break;
		case 'W':
			EBCDIC_str[i] = 230;
			break;
		case 'X':
			EBCDIC_str[i] = 231;
			break;
		case 'Y':
			EBCDIC_str[i] = 232;
			break;
		case 'Z':
			EBCDIC_str[i] = 233;
			break;
		case '0':
			EBCDIC_str[i] = 240;
			break;
		case '1':
			EBCDIC_str[i] = 241;
			break;
		case '2':
			EBCDIC_str[i] = 242;
			break;
		case '3':
			EBCDIC_str[i] = 243;
			break;
		case '4':
			EBCDIC_str[i] = 244;
			break;
		case '5':
			EBCDIC_str[i] = 245;
			break;
		case '6':
			EBCDIC_str[i] = 246;
			break;
		case '7':
			EBCDIC_str[i] = 247;
			break;
		case '8':
			EBCDIC_str[i] = 248;
			break;
		case '9':
			EBCDIC_str[i] = 249;
			break;
		default:
			EBCDIC_str[i] = 64;  /* space */
			break;
		}
	}
}

