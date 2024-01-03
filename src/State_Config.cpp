
#include "State_Config.hpp"


namespace cuttlefish
{

uint8_t State_Config::f_th;


void State_Config::set_edge_threshold(const uint8_t f_th)
{
    State_Config::f_th = f_th;
}

}
