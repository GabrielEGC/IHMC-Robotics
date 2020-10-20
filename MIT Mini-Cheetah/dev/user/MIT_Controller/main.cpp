/*!
 * @file main.cpp
 * @brief Main Function for the WBC Controller
 *
 * The main function parses command line arguments and starts the appropriate
 * driver.
 */

#include <main_helper.h>
#include "MIT_Controller.hpp"

int main(int argc, char** argv) {
  main_helper(argc, argv, new MIT_Controller());
  return 0;
}
