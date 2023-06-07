#include "s21_decimal.h"
#include <stdio.h>

int s21_add(s21_decimal value_1, s21_decimal value_2, s21_decimal *result) {
    int status = 0;
    int sign1, sign2;
    int check = 0;

    sign1 = getSign(value_1);
    sign2 = getSign(value_2);

    setSign(&value_1, 0);
    setSign(&value_2, 0);

    initDecimal(result);
    convertorScale(&value_1, &value_2);

    while (check == 0) {
        if ((sign1 == 0 && sign2 == 0) || (sign1 == 1 && sign2 == 1)) {
            status = sumBits(value_1, value_2, result);
            if (sign1 == 1 && status != 0) {
                status = 2;
            }
            setSign(result, sign1);
        } else {
            if (s21_is_less(value_1, value_2)) {
                status = subBits(value_2, value_1, result);
                setSign(result, sign2);
            } else {
                status = subBits(value_1, value_2, result);
                setSign(result, sign1);
            }
        }
        if (status == 0) {
            check = 1;
        } else if (getScale(value_1) == 0) {
            check = 1;
        } else {
            shiftScaleRigth(&value_1, 1, 0);
            shiftScaleRigth(&value_2, 1, 0);
        }
    }
    setScale(result, getScale(value_1));
    return status;
}

int s21_sub(s21_decimal value_1, s21_decimal value_2, s21_decimal *result) {
    int status = 0;
    int sign1, sign2;
    int check = 0;

    sign1 = getSign(value_1);
    sign2 = getSign(value_2);

    setSign(&value_1, 0);
    setSign(&value_2, 0);

    initDecimal(result);
    convertorScale(&value_1, &value_2);

    while (check == 0) {
        if ((sign1 == 0 && sign2 == 0)) {
            if (s21_is_less_or_equal(value_1, value_2)) {
                status = subBits(value_2, value_1, result);
                setSign(result, 1);
            } else {
                status = subBits(value_1, value_2, result);
                setSign(result, 0);
            }
        } else if ((sign1 == 1 && sign2 == 1)) {
            if (s21_is_less_or_equal(value_1, value_2)) {
                status = subBits(value_2, value_1, result);
                setSign(result, 0);
            } else {
                status = subBits(value_1, value_2, result);
                setSign(result, 1);
            }
        } else if ((sign1 == 0 && sign2 == 1)) {
            status = sumBits(value_1, value_2, result);
            setSign(result, 0);
        } else if ((sign1 == 1 && sign2 == 0)) {
            status = sumBits(value_1, value_2, result);
            setSign(result, 1);
            if (status != 0) {
                status = 2;
            }
        }
        if (status == 0) {
            check = 1;
        } else if (getScale(value_1) == 0) {
            check = 1;
        } else {
            shiftScaleRigth(&value_1, 1, 0);
            shiftScaleRigth(&value_2, 1, 0);
        }
    }
    setScale(result, getScale(value_1));
    return status;
}

int s21_mul(s21_decimal value_1, s21_decimal value_2, s21_decimal *result) {
    s21_decimal tmp;
    int status = 0;

    initDecimal(result);
    for (int i = 0; i <= positionFirstNum(value_2); i++) {
        if (getBit(value_2, i) == 1) {
            initDecimal(&tmp);
            tmp = value_1;
            if (shiftDecimalLeft(&tmp, i) == 0) {
                if (sumBits(*result, tmp, result)) {
                    status = choiceForShiftDecimalRight(result, &value_1,
                                                        &value_2, &i);
                }
            } else {
                status =
                    choiceForShiftDecimalRight(result, &value_1, &value_2, &i);
            }
        }
    }

    setSignForResult(result, value_1, value_2);
    if (status == 0) {
        int Scale1, Scale2, Scale;
        Scale1 = getScale(value_1);
        Scale2 = getScale(value_2);
        Scale = Scale1 + Scale2;
        if (Scale > 28) {
            status = 2;
        } else {
            setScale(result, Scale);
        }
    }
    return status;
}

int s21_div(s21_decimal value_1, s21_decimal value_2, s21_decimal *result) {
    int status = 0;
    int sign1 = getSign(value_1);
    int sign2 = getSign(value_2);

    if (checkNull(value_2) == 0) {
        status = 3;
    } else if (checkNull(value_1) == 0) {
        initDecimal(result);
    } else {
        convertorScale(&value_1, &value_2);
        if (checkNull(value_1) == 0) {
            status = 2;
        } else if (checkNull(value_2) == 0) {
            status = 1;
            if (sign1 != sign2) {
                status = 2;
            }
        } else {
            divOfNumbers(&value_1, &value_2, result);
            if (checkNull(*result) == 0) {
                status = 2;
            }
        }
        if (sign1 != sign2) {
            setSign(result, 1);
        }
    }
    return status;
}

int s21_mod(s21_decimal value_1, s21_decimal value_2, s21_decimal *result) {
    int status;
    int sign1 = getSign(value_1);
    int sign2 = getSign(value_2);
    s21_decimal incomplete_quotient;

    setSign(&value_1, 0);
    setSign(&value_2, 0);
    status = s21_div(value_1, value_2, result);
    if (status == 1) {
        result->bits[0] = 0;
        result->bits[1] = 0;
        result->bits[2] = 0;
        result->bits[3] = 0;
        status = 0;
    } else if (status == 0) {
        if (s21_is_less(value_1, value_2)) {
            *result = value_1;
        } else {
            convertorScale(&value_1, &value_2);
            s21_truncate(*result, &incomplete_quotient);
            s21_mul(value_2, incomplete_quotient, &value_2);
            s21_sub(value_1, value_2, result);
        }
    } else if (status == 2) {
        *result = value_1;
        status = 0;
    }
    if (sign1 == 1 && sign2 == 1) {
        setSign(result, 1);
    }
    return status;
}

int s21_is_less(s21_decimal value_1, s21_decimal value_2) {
    int status;
    int sign1 = getSign(value_1);
    int sign2 = getSign(value_2);
    int check = 0;

    convertorScale(&value_1, &value_2);
    if (sign1 == 1 && sign2 == 0) {
        status = 1;
    } else if (sign1 == 0 && sign2 == 1) {
        status = 0;
    } else {
        for (int i = 95; i >= 0 && check == 0; i--) {
            if (getBit(value_1, i) < getBit(value_2, i)) {
                check = 1;
            } else if (getBit(value_1, i) > getBit(value_2, i)) {
                check = 2;
            }
        }
        if ((check == 1 && sign1 == 0 && sign2 == 0) ||
            (check == 2 && sign1 == 1 && sign2 == 1)) {
            status = 1;
        } else {
            status = 0;
        }
    }
    if (value_1.bits[0] == 0 && value_2.bits[0] == 0 && value_1.bits[1] == 0 &&
        value_2.bits[1] == 0 && value_1.bits[2] == 0 && value_2.bits[2] == 0) {
        status = 0;
    }
    return status;
}

int s21_is_less_or_equal(s21_decimal value_1, s21_decimal value_2) {
    return s21_is_less(value_1, value_2) || s21_is_equal(value_1, value_2);
}

int s21_is_greater(s21_decimal value_1, s21_decimal value_2) {
    return !s21_is_less_or_equal(value_1, value_2);
}

int s21_is_greater_or_equal(s21_decimal value_1, s21_decimal value_2) {
    return !s21_is_less(value_1, value_2);
}

int s21_is_equal(s21_decimal value_1, s21_decimal value_2) {
    int status = 1;
    int sign1 = getSign(value_1);
    int sign2 = getSign(value_2);

    convertorScale(&value_1, &value_2);
    if (sign1 == sign2) {
        for (int i = 95; i >= 0 && status == 1; i--) {
            if (getBit(value_1, i) != getBit(value_2, i)) {
                status = 0;
            }
        }
    } else {
        status = 0;
    }
    if (value_1.bits[0] == 0 && value_2.bits[0] == 0 && value_1.bits[1] == 0 &&
        value_2.bits[1] == 0 && value_1.bits[2] == 0 && value_2.bits[2] == 0) {
        status = 1;
    }
    return status;
}

int s21_is_not_equal(s21_decimal value_1, s21_decimal value_2) {
    return !s21_is_equal(value_1, value_2);
}

int s21_from_int_to_decimal(int src, s21_decimal *dst) {
    int status = 0;
    if (dst) {
        initDecimal(dst);
        if (src >= 0) {
            dst->bits[0] = src;
        } else {
            setSign(dst, 1);
            dst->bits[0] = ~src + 1;
        }
    } else {
        status = 1;
    }
    return status;
}

int s21_from_float_to_decimal(float src, s21_decimal *dst) {
    int status = 0;
    if (src == 1 / 0.0 && src == 0 / 0.0 && src == -1 / 0.0 && src == -0 / 0.0)
        status = 1;
    initDecimal(dst);
    if (src != 0) {
        initDecimal(dst);
        if (dst && !status) {
            int is_negative = 0;
            if (src < 0) {
                is_negative = 1;
                src *= -1;
            }
            double dbl = src;
            char float_bin_buff[32] = "";
            getBinaryFromFloat((float)dbl, float_bin_buff);
            int exp = getScaleFloat(float_bin_buff);
            int scale = 0;
            while (scale < 28 && (int)dbl / (int)pow(2, 21) == 0) {
                dbl *= 10;
                scale++;
            }
            dbl = round(dbl);
            if (scale <= 28 && (exp > -94 && exp < 96)) {
                while (fmod(dbl, 10) == 0 && scale > 0) {
                    dbl = dbl / 10;
                    scale--;
                }
                getBinaryFromFloat((float)dbl, float_bin_buff);
                exp = getScaleFloat(float_bin_buff);
                setBit(dst, exp, 1);
                writeMantissaToDecimal(dst, float_bin_buff, exp);
                setScale(dst, scale);
                setSign(dst, is_negative);
            } else {
                status = 1;
            }
        }
    }
    return status;
}

int s21_from_decimal_to_int(s21_decimal src, int *dst) {
    int status = 0;
    if (!dst) {
        status = 1;
    } else {
        s21_decimal truncated;
        status = s21_truncate(src, &truncated);
        int highBit = positionFirstNum(truncated);
        if (highBit > 31) {
            status = 1;
        } else {
            int base = 1;
            *dst = 0;
            for (int i = 0; i <= highBit; i++) {
                *dst += getBit(truncated, i) * base;
                base = base * 2;
            }
            if (getBit(src, 127)) *dst *= -1;
        }
    }
    return status;
}

int s21_from_decimal_to_float(s21_decimal src, float *dst) {
    int ret = 0;
    if (!dst || getScale(src) > 28) {
        ret = 1;
    } else {
        *dst = 0.0;
        int scale = getScale(src);
        int last_bit = getBit(src, 127);
        int sign = getBit(src, 127);
        if (sign) setBit(&src, 127, 1);
        setScale(&src, scale);
        double tmp = *dst;
        if (last_bit == 1) setBit(&src, 127, 1);
        unsigned long base = 1;
        for (int i = 0; i < 96; i++) {
            tmp += getBit(src, i) * base;
            base = base * 2;
        }
        while (scale != 0) {
            tmp = tmp / 10;
            scale--;
        }
        if (getBit(src, 127)) tmp *= -1;
        *dst = tmp;
    }
    return ret;
}

int s21_floor(s21_decimal value, s21_decimal *result) {
    int Scale = getScale(value);
    int Sign = getSign(value);
    int check = decimalCheck(value);

    shiftScaleRigth(&value, Scale, 1);
    if (check == 1 && Sign == 1) {
        sumWithNumOne(&value);
    }
    *result = value;
    setSign(result, Sign);
    return 0;
}

int s21_round(s21_decimal value, s21_decimal *result) {
    int Scale = getScale(value);
    int Sign = getSign(value);

    shiftScaleRigth(&value, Scale, 0);
    *result = value;
    setSign(result, Sign);
    return 0;
}

int s21_truncate(s21_decimal value, s21_decimal *result) {
    int Scale = getScale(value);
    int Sign = getSign(value);

    shiftScaleRigth(&value, Scale, 1);
    *result = value;
    setSign(result, Sign);
    return 0;
}

int s21_negate(s21_decimal value, s21_decimal *result) {
    int Sign = getSign(value);

    *result = value;
    if (Sign == 0) {
        setSign(result, 1);
    } else {
        setSign(result, 0);
    }
    return 0;
}

void setBit(s21_decimal *d, int i, int v) {
    unsigned int mask = 1u << (i % 32);
    if (v == 1) {
        d->bits[i / 32] |= mask;
    }
    if (v == 0) {
        d->bits[i / 32] &= ~mask;
    }
}

int getBit(s21_decimal d, int i) {
    unsigned int mask = 1u << (i % 32);
    return ((d.bits[i / 32] & mask) != 0);
}

int getSign(s21_decimal d) { return d.bits[3] >> 31; }

void setSign(s21_decimal *d, int i) {
    unsigned int mask = 1u << 31;
    if (i == 1) {
        d->bits[3] |= mask;
    }
    if (i == 0) {
        d->bits[3] &= ~mask;
    }
}

int getScale(s21_decimal d) {
    setSign(&d, 0);
    d.bits[3] >>= 16;
    return d.bits[3];
}

void setScale(s21_decimal *d, int i) {
    if (i >= 0 && i <= 28) {
        int cleanMask = ~(0xFF << 16);
        d->bits[3] &= cleanMask;
        int mask = i << 16;
        d->bits[3] |= mask;
    }
}

int shiftScaleLeft(s21_decimal *d, int scale) {
    int status = 0;
    int origin_sign = getScale(*d);
    s21_decimal tmp_result = *d;
    setScale(&tmp_result, 0);
    for (int i = 0; i < scale; i++) {
        s21_decimal tmp1 = tmp_result;
        s21_decimal tmp2 = tmp_result;
        if (shiftDecimalLeft(&tmp1, 3) == 1) {
            status = 1;
            break;
        }
        shiftDecimalLeft(&tmp2, 1);
        if (sumBits(tmp1, tmp2, &tmp_result) == 1) {
            status = 1;
            break;
        }
        if (status == 0) {
            origin_sign++;
            initDecimal(d);
            *d = tmp_result;
            setScale(d, origin_sign);
        }
    }
    return status;
}

int shiftDecimalLeft(s21_decimal *d, int shift) {
    int status = 0;
    for (int i = 0; i < shift; i++) {
        int bits0 = getBit(*d, 31);
        int bits1 = getBit(*d, 63);
        int bits2 = getBit(*d, 95);
        d->bits[0] <<= 1;
        d->bits[1] <<= 1;
        d->bits[2] <<= 1;
        if (bits0 == 1) {
            setBit(d, 32, 1);
        }
        if (bits1 == 1) {
            setBit(d, 64, 1);
        }
        if (bits2 == 1) {
            status = 1;
            break;
        }
    }
    return status;
}

int shiftScaleRigth(s21_decimal *d, int shiftScale, int status) {
    unsigned int minuend;
    s21_decimal resultOfDevide;
    int checkForMinuend = 0;

    int Scale = getScale(*d) - shiftScale;
    for (; shiftScale != 0; shiftScale--) {
        initDecimal(&resultOfDevide);
        int begin = positionFirstNum(*d);

        minuend = 0;
        initMinuend(&minuend, *d, begin);
        int i = 95;
        int check = 0;

        for (int j = begin - 3; j >= 0 && check == 0; j--, i--) {
            if (minuend >= 10u) {
                setBit(&resultOfDevide, i, 1);
                minuend -= 10u;
                minuend <<= 1;
                setBitInt(&minuend, 0, getBit(*d, j - 1));
            } else {
                setBit(&resultOfDevide, i, 0);
                minuend <<= 1;
                if (j != 0) {
                    setBitInt(&minuend, 0, getBit(*d, j - 1));
                }
            }
            if (minuend < 10u && j == 0) {
                check = 1;
            }
        }
        minuend >>= 1;
        initDecimal(d);
        int tmp = 95;
        for (int j = 94 - i; j >= 0; j--, tmp--) {
            setBit(d, j, getBit(resultOfDevide, tmp));
        }
        if (status == 0) {
            if (minuend != 0 && shiftScale > 1) {
                checkForMinuend++;
            }
            if (minuend >= 5 && shiftScale == 1 && checkForMinuend != 0) {
                sumWithNumOne(d);
            }
            if (minuend == 5 && shiftScale == 1 && checkForMinuend == 0) {
                int even = evenOrNot(*d);
                if (even == 1) {
                    sumWithNumOne(d);
                }
            }
        }
    }
    setScale(d, Scale);
    return 0;
}

void initMinuend(unsigned int *minuend, s21_decimal d, int begin) {
    int i = 3;
    for (int j = 0; i >= 0; j++, i--) {
        setBitInt(minuend, j, getBit(d, begin - i));
    }
}

void sumWithNumOne(s21_decimal *d) {
    int changer = 0;
    s21_decimal value_1 = *d;
    s21_decimal value_2 = {{1, 0, 0, 0}};

    initDecimal(d);
    shiftScaleLeft(&value_2, getScale(value_1));
    for (int i = 0; i < 96; i++) {
        int num_1 = getBit(value_1, i);
        int num_2 = getBit(value_2, i);
        if ((num_1 == 0 && num_2 == 0) || (num_1 == 1 && num_2 == 1))
            setBit(d, i, 0 + changer);
        if ((num_1 == 0 && num_2 == 1) || (num_1 == 1 && num_2 == 0))
            setBit(d, i, 1 - changer);
        if ((num_1 + num_2 + changer) > 1)
            changer = 1;
        else
            changer = 0;
    }
}

void setBitInt(unsigned int *d, int i, int v) {
    unsigned int mask = 1u << i;
    if (v == 1) *d |= mask;
    if (v == 0) *d &= ~mask;
}

int convertorScale(s21_decimal *value_1, s21_decimal *value_2) {
    int Scale1 = getScale(*value_1);
    int Scale2 = getScale(*value_2);
    int difference;
    if (Scale1 > Scale2) {
        difference = Scale1 - Scale2;
        if (shiftScaleLeft(value_2, difference) == 1) {
            Scale2 = getScale(*value_2);
            difference = Scale1 - Scale2;
            shiftScaleRigth(value_1, difference, 0);
        }
    } else if (Scale1 < Scale2) {
        difference = Scale2 - Scale1;
        if (shiftScaleLeft(value_1, difference) == 1) {
            Scale1 = getScale(*value_1);
            difference = Scale2 - Scale1;
            shiftScaleRigth(value_2, difference, 0);
        }
    }
    return 0;
}

int sumBits(s21_decimal value_1, s21_decimal value_2, s21_decimal *result) {
    int changer = 0;
    int status = 0;

    for (int i = 0; i < 96; i++) {
        int num_1 = getBit(value_1, i);
        int num_2 = getBit(value_2, i);
        if ((num_1 == 0 && num_2 == 0) || (num_1 == 1 && num_2 == 1)) {
            setBit(result, i, 0 + changer);
        } else if ((num_1 == 0 && num_2 == 1) || (num_1 == 1 && num_2 == 0)) {
            setBit(result, i, 1 - changer);
        }
        if ((num_1 + num_2 + changer) > 1) {
            changer = 1;
        } else {
            changer = 0;
        }
        if (i == 95 && changer == 1) {
            status = 1;
        }
    }
    return status;
}

int subBits(s21_decimal value_1, s21_decimal value_2, s21_decimal *result) {
    int big_value = 0;

    for (int i = 0; i < 96; i++) {
        int num_1 = getBit(value_1, i);
        int num_2 = getBit(value_2, i);
        if ((num_1 == 0 && num_2 == 0) || (num_1 == 1 && num_2 == 1)) {
            setBit(result, i, 0 + big_value);
        } else if ((num_1 == 0 && num_2 == 1) || (num_1 == 1 && num_2 == 0)) {
            setBit(result, i, 1 - big_value);
        }
        if ((num_1 - num_2 - big_value) < 0) {
            big_value = 1;
        } else {
            big_value = 0;
        }
    }
    return 0;
}

int positionFirstNum(s21_decimal d) {
    int i = 95;
    for (; getBit(d, i) != 1 && i >= 0; i--) {
    }
    return i;
}

void initDecimal(s21_decimal *d) {
    if (d) {
        d->bits[0] = 0;
        d->bits[1] = 0;
        d->bits[2] = 0;
        d->bits[3] = 0;
    }
}

int getScaleFloat(const char *src) {
    int result = 0, base = 1;
    for (int i = 8; i > 0; i--) {
        result += src[i] * base;
        base = base * 2;
    }
    return (result - 127);
}

void getBinaryFromFloat(float src, char *float_bin_buff) {
    union float_value float_number_bits;
    float_number_bits.float_view = src;
    for (unsigned int mask = 0x80000000; mask; mask >>= 1) {
        *float_bin_buff = !!(float_number_bits.int_view & mask);
        float_bin_buff++;
    }
}

void writeMantissaToDecimal(s21_decimal *d, char *float_bin_buff, int exp) {
    for (int i = exp - 1, j = 9; i >= 0 && j < 32; i--, j++) {
        if (float_bin_buff[j]) setBit(d, i, 1);
    }
}

int choiceForShiftDecimalRight(s21_decimal *result, s21_decimal *value_1,
                               s21_decimal *value_2, int *i) {
    int Scale1, Scale2, Sign1, Sign2;
    int status = 0;

    Scale1 = getScale(*value_1);
    Scale2 = getScale(*value_2);
    Sign1 = getSign(*value_1);
    Sign2 = getSign(*value_2);
    if (Scale1 != 0 || Scale2 != 0) {
        if (Scale1 >= Scale2) {
            shiftScaleRigth(value_1, 1, 0);
            setSign(value_1, Sign1);
        } else {
            shiftScaleRigth(value_2, 1, 0);
            setSign(value_2, Sign2);
        }
        initDecimal(result);
        *i = -1;
    } else {
        if (Sign1 == Sign2) {
            status = 1;
        } else {
            status = 2;
        }
        *i = 100;
    }
    return status;
}

void setSignForResult(s21_decimal *result, s21_decimal value_1,
                      s21_decimal value_2) {
    int sign1 = getSign(value_1);
    int sign2 = getSign(value_2);

    if ((sign1 && !sign2) || (!sign1 && sign2)) {
        setSign(result, 1);
    }
}

int decimalCheck(s21_decimal d) {
    int status = 1;
    int Scale = getScale(d);
    s21_decimal tmp;
    s21_decimal result;

    initDecimal(&tmp);
    initDecimal(&result);
    if (Scale != 0) {
        s21_truncate(d, &tmp);
        s21_sub(d, tmp, &result);
        if (result.bits[0] == 0 && result.bits[1] == 0 && result.bits[2] == 0) {
            status = 0;
        }
    } else {
        status = 0;
    }
    return status;
}

int evenOrNot(s21_decimal d) {
    unsigned int minuend;
    int begin = positionFirstNum(d);
    int result = 1;

    minuend = 0;
    if (begin >= 3) {
        initMinuend(&minuend, d, begin);
        int i = 95;
        int check = 0;

        for (int j = begin - 3; j >= 0 && check == 0; j--, i--) {
            if (minuend >= 10u) {
                minuend -= 10u;
                minuend <<= 1;
                setBitInt(&minuend, 0, getBit(d, j - 1));
            } else {
                minuend <<= 1;
                if (j != 0) {
                    setBitInt(&minuend, 0, getBit(d, j - 1));
                }
            }
            if (minuend < 10u && j == 0) {
                check = 1;
            }
        }
        minuend >>= 1;
    } else {
        minuend = d.bits[0];
    }
    if (minuend % 2 == 0 || minuend == 0) {
        result = 0;
    }
    return result;
}

void divOfNumbers(s21_decimal *value_1, s21_decimal *value_2,
                  s21_decimal *result) {
    s21_decimal minuend, resultOfDevide;
    int i, j, status_of_mul;
    int firstCycle = 0;
    s21_decimal ten = {{10, 0, 0, 0}};
    int wholeNumber = 0;

    do {
        status_of_mul = 0;
        initDecimal(&minuend);
        initDecimal(&resultOfDevide);
        nullSignAndScale(value_1, value_2);
        if (s21_is_less(*value_1, *value_2)) {
            minuend = *value_1;
            i = -1;
            if (firstCycle == 1) {
                wholeNumber--;
            }
        } else {
            i = initDecimalMinuend(&minuend, *value_1, *value_2);
        }

        j = 95;
        int check = 0;
        for (; j >= 0 && check == 0; i--, j--) {
            if (s21_is_greater_or_equal(minuend, *value_2)) {
                setBit(&resultOfDevide, j, 1);
                subBits(minuend, *value_2, &minuend);
                if (i >= 0) {
                    shiftDecimalLeft(&minuend, 1);
                    setBit(&minuend, 0, getBit(*value_1, i));
                }
            } else {
                setBit(&resultOfDevide, j, 0);
                if (i >= 0) {
                    shiftDecimalLeft(&minuend, 1);
                    setBit(&minuend, 0, getBit(*value_1, i));
                }
            }
            if (s21_is_less(minuend, *value_2) && i == -1 && firstCycle == 0) {
                check = 1;
                if (checkNull(resultOfDevide) != 0) {
                    initDecimal(result);
                    int tmp = 95;
                    for (int k = 95 - j; k >= 0; k--, tmp--) {
                        setBit(result, k, getBit(resultOfDevide, tmp));
                    }
                    wholeNumber = sizeOfNumber(*result);
                    initDecimal(result);
                }
                firstCycle = 1;
            }
            if (s21_is_less(minuend, *value_2) && i == -1 && firstCycle == 1) {
                check = 1;
            }
        }
        status_of_mul = s21_mul(*value_1, ten, value_1);
    } while (checkNull(minuend) != 0 && status_of_mul == 0);
    initDecimal(result);
    int tmp = 95;
    for (i = 94 - j; i >= 0; i--, tmp--) {
        setBit(result, i, getBit(resultOfDevide, tmp));
    }
    int fractionNumber = sizeOfNumber(*result);
    int Scale = fractionNumber - wholeNumber;
    setScale(result, Scale);
}

int initDecimalMinuend(s21_decimal *minuend, s21_decimal value_1,
                       s21_decimal value_2) {
    int i = positionFirstNum(value_1);
    s21_decimal tmp_minuend = *minuend;
    initDecimal(&tmp_minuend);
    while (s21_is_less(tmp_minuend, value_2)) {
        shiftDecimalLeft(&tmp_minuend, 1);
        setBit(&tmp_minuend, 0, getBit(value_1, i));
        i--;
    }
    *minuend = tmp_minuend;
    return i;
}

void nullSignAndScale(s21_decimal *d1, s21_decimal *d2) {
    setScale(d1, 0);
    setScale(d2, 0);
    setSign(d1, 0);
    setSign(d2, 0);
}

int checkNull(s21_decimal d) {
    int status = 1;
    if (d.bits[0] == 0 && d.bits[1] == 0 && d.bits[2] == 0) {
        status = 0;
    }
    return status;
}

int sizeOfNumber(s21_decimal d) {
    s21_decimal resultOfDevide;
    unsigned int minuend;
    int size = 0;

    while (checkNull(d) != 0) {
        initDecimal(&resultOfDevide);
        int begin = positionFirstNum(d);

        minuend = 0;
        initMinuend(&minuend, d, begin);
        int i = 95;
        int check = 0;

        for (int j = begin - 3; j >= 0 && check == 0; j--, i--) {
            if (minuend >= 10u) {
                setBit(&resultOfDevide, i, 1);
                minuend -= 10u;
                minuend <<= 1;
                setBitInt(&minuend, 0, getBit(d, j - 1));
            } else {
                setBit(&resultOfDevide, i, 0);
                minuend <<= 1;
                if (j != 0) {
                    setBitInt(&minuend, 0, getBit(d, j - 1));
                }
            }
            if (minuend < 10u && j == 0) {
                check = 1;
            }
        }
        initDecimal(&d);
        int tmp = 95;
        for (int j = 94 - i; j >= 0; j--, tmp--) {
            setBit(&d, j, getBit(resultOfDevide, tmp));
        }
        size++;
    }
    return size;
}
