
const char* rst = "\033[m";
const char* bld = "\033[1m";
const char* red = "\033[31m";
const char* grn = "\033[32m";
const char* yel = "\033[33m";
const char* blu = "\033[34m";
const char* mag = "\033[35m";
const char* cya = "\033[36m";


const DeviceEntry my_devices[] = {
    DeviceEntry(0x11111111, "leak sensor laundry",      red),
    DeviceEntry(0x22222222, "leak sensor refrigerator", grn),
    DeviceEntry(0x33333333, "leak sensor crawlspace",   yel),

    DeviceEntry(0x44444444, "temp sensor office",       blu),
    DeviceEntry(0x55555555, "temp sensor refrigerator", mag),
    DeviceEntry(0x66666666, "temp sensor freezer",      cya),
};

const uint32 my_devices_len = (sizeof(my_devices) / sizeof(my_devices[0]));
