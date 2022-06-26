#ifndef SRC_S21_DECIMAL_H_
#define SRC_S21_DECIMAL_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    unsigned int bits[4];
} s21_decimal;

union float_value {
    unsigned int int_view;
    float float_view;
};

int s21_add(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int s21_sub(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int s21_mul(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int s21_div(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int s21_mod(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);

int s21_is_less(s21_decimal value_1, s21_decimal value_2);
int s21_is_less_or_equal(s21_decimal value_1, s21_decimal value_2);
int s21_is_greater(s21_decimal value_1, s21_decimal value_2);
int s21_is_greater_or_equal(s21_decimal value_1, s21_decimal value_2);
int s21_is_equal(s21_decimal value_1, s21_decimal value_2);
int s21_is_not_equal(s21_decimal value_1, s21_decimal value_2);

int s21_from_int_to_decimal(int src, s21_decimal *dst);
int s21_from_float_to_decimal(float src, s21_decimal *dst);
int s21_from_decimal_to_int(s21_decimal src, int *dst);
int s21_from_decimal_to_float(s21_decimal src, float *dst);

int s21_floor(s21_decimal value, s21_decimal *result);
int s21_round(s21_decimal value, s21_decimal *result);
int s21_truncate(s21_decimal value, s21_decimal *result);
int s21_negate(s21_decimal value, s21_decimal *result);

void setBit(s21_decimal *d, int i, int v);
int getBit(s21_decimal d, int i);
int getSign(s21_decimal d);
void setSign(s21_decimal *d, int i);
int getScale(s21_decimal d);
void setScale(s21_decimal *d, int i);
int isNegative(s21_decimal *value_1, s21_decimal *value_2);
int convertorScale(s21_decimal *value_1, s21_decimal *value_2);
int shiftDecimalLeft(s21_decimal *d, int shift);
int shiftScaleLeft(s21_decimal *d, int scale);
int sumBits(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int subBits(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int shiftScaleRigth(s21_decimal *d, int shiftScale, int status);
void initMinuend(unsigned int *minuend, s21_decimal d, int begin);
void sumWithNumOne(s21_decimal *d);
void setBitInt(unsigned int *d, int i, int v);
int positionFirstNum(s21_decimal d);
void initDecimal(s21_decimal *d);
int getScaleFloat(const char *src);
void getBinaryFromFloat(float src, char *float_bin_buff);
void writeMantissaToDecimal(s21_decimal *d, char *float_bin_buff, int exp);
int choiceForShiftDecimalRight(s21_decimal *result, s21_decimal *value_1,
                               s21_decimal *value_2, int *i);
void setSignForResult(s21_decimal *result, s21_decimal value_1,
                      s21_decimal value_2);
int decimalCheck(s21_decimal d);
int evenOrNot(s21_decimal d);
int initDecimalMinuend(s21_decimal* minuend, s21_decimal value_1, s21_decimal value_2);
void nullSignAndScale(s21_decimal* d1, s21_decimal* d2);
int checkNull(s21_decimal d);
int sizeOfNumber(s21_decimal d);
void divOfNumbers(s21_decimal* value_1, s21_decimal* value_2, s21_decimal* result);

#endif  // SRC_S21_DECIMAL_H_
