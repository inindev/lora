

class DeviceEntry {
private:
    const uint32_t addr;
    const char* name;
    const char* color;
    mutable uint32_t last_millis;
public:
    DeviceEntry(const uint32_t addr, const char* name, const char* color="") : addr(addr), name(name), color(color), last_millis(0) { }
    const uint32_t get_addr() const { return this->addr; }
    const char* get_name() const { return this->name; }
    const char* get_color() const { return this->color; }
    void set_last_seen() const { this->last_millis = millis(); }
    const uint32_t get_elapsed_sec() const {
    	if(0 == this->last_millis) return 0;
    	return (millis() - this->last_millis) / 1000;
    }
};

#include "my_devices.h"

const DeviceEntry* get_device_entry(const uint32_t dev_addr) {
    const uint32_t len = my_devices_len;
    for(uint32_t i=0; i<len; i++) {
        if(dev_addr == my_devices[i].get_addr()) {
            return &my_devices[i];
        }
    }
    return NULL;
}

