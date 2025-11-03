#include <stdio.h>
#include <config/ConfigParser.hpp>

int main() {
  ConfigParser().parse("../config.yml", "soda1");
}