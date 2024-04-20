
// Copyright (C) 2024, John Clark <inindev@gmail.com>

#include <stdio.h>
#include <stdint.h>
#include "device_entry.h"


typedef enum {
    mtype_join_request          = 0,
    mtype_join_accept           = 1,
    mtype_unconfirmed_data_up   = 2,
    mtype_unconfirmed_data_down = 3,
    mtype_confirmed_data_up     = 4,
    mtype_confirmed_data_down   = 5,
    mtype_rejoin_request        = 6,
    mtype_proprietary           = 7,
} msg_type_t;

typedef enum {
    mver_lorawan_v1 = 0,
} msg_ver_t;


uint8_t get_mhdr(const uint8_t* payload) {
    return (payload[0]);
}
msg_type_t get_msg_type(const uint8_t* payload) {
    return (msg_type_t)(payload[0] >> 5);
}
msg_ver_t get_msg_ver(const uint8_t* payload) {
    return (msg_ver_t)(payload[0] & 0x03);
}

const char* get_msg_dir(const msg_type_t msg_type) {
    if(msg_type > mtype_rejoin_request) return "";
    if(msg_type % 2 == 0) return "up";
    return "down";
}

const char* get_msg_type_str(const msg_type_t msg_type) {
    switch(msg_type) {
    case mtype_join_request:
        return "join request";
    case mtype_join_accept:
        return "join accept";
    case mtype_unconfirmed_data_up:
        return "unconfirmed data up";
    case mtype_unconfirmed_data_down:
        return "unconfirmed data down";
    case mtype_confirmed_data_up:
        return "confirmed data up";
    case mtype_confirmed_data_down:
        return "confirmed data down";
    case mtype_rejoin_request:
        return "rejoin request";
    case mtype_proprietary:
        return "proprietary";
    default:
        return "unknown mhdr type";
    }
}

const char* get_msg_ver_str(const msg_ver_t msg_ver) {
    switch(msg_ver) {
    case mver_lorawan_v1:
        return "lorawan v1";
    default:
        return "unknown mhdr version";
    }
}

// last 4 bytes of payload
const uint8_t* get_mic(const uint8_t* payload, const uint8_t payload_len) {
    return (&payload[payload_len - 4]);
}

void print_hex(const uint8_t* buf, const uint8_t len) {
    for(uint8_t i=0; i<len; i++) {
        printf("%02hhx", buf[i]);
    }
    printf("\r\n");
}

void print_mac_payload(const uint8_t* payload, const uint8_t payload_len) {
    print_hex(&payload[1], payload_len - 5);
}

// first 7 bytes of mac payload
const uint8_t* get_fhdr(const uint8_t* payload) {
    return (&payload[1]);
}

const uint8_t get_fport(const uint8_t* payload) {
    return (payload[8]);
}

void print_frm_payload(const uint8_t* payload, const uint8_t payload_len) {
    print_hex(&payload[9], payload_len - 13);
}

uint32_t get_dev_addr(const uint8_t* payload) {
    return(payload[4] << 24 | payload[3] << 16 | payload[2] << 8 | payload[1]);
}

uint8_t get_fctrl(const uint8_t* payload) {
    return (payload[5]);
}

uint16_t get_fcnt(const uint8_t* payload) {
    return (payload[7] << 8 | payload[6]);
}

// fopts[0..15] ??

/////////////////////////////////
void print_mini_report(const uint8_t* payload, const uint8_t payload_len, const bool colorize=false) {
    const uint32_t dev_addr = get_dev_addr(payload);
    const DeviceEntry* device_entry = get_device_entry(dev_addr);

    const char* dev_name  = (NULL == device_entry) ? "" : device_entry->get_name();
    const char* dev_color = (NULL == device_entry) ? "" : device_entry->get_color();
    const uint32_t el_sec = (NULL == device_entry) ?  0 : device_entry->get_elapsed_sec();
    if(NULL != device_entry) device_entry->set_last_seen();

    if(colorize) printf(dev_color);
    printf("    device addr: %s (%08x)\r\n", dev_name, dev_addr);
    printf("   message type: %s (0x%02hhx)\r\n", get_msg_type_str(get_msg_type(payload)), get_msg_type(payload));
    printf("      direction: %s\r\n", get_msg_dir(get_msg_type(payload)));
    printf("           fcnt: %d\r\n", get_fcnt(payload));
    printf("    phy payload: "); print_hex(payload, payload_len);
    printf("    frm payload:                   "); print_frm_payload(payload, payload_len);
    if(el_sec > 0) printf("      last seen: %ds ago\r\n", el_sec);
    if(colorize) printf(rst);
}

void print_report(const uint8_t* payload, const uint8_t payload_len) {
    printf("---------\r\n");
    printf("phy payload: ");
    print_hex(payload, payload_len);

    printf("---------\r\n");
    printf("mhdr: %02hhx\n", get_mhdr(payload));

    printf("mac payload: ");
    print_mac_payload(payload, payload_len);

    const uint8_t* mic = get_mic(payload, payload_len);
    printf("mic: ");
    print_hex(mic, 4); // mic is 4 bytes

    printf("---------\r\n");
    const uint8_t* fhdr = get_fhdr(payload);
    printf("fhdr: ");
    print_hex(fhdr, 7); // fhdr is 7 bytes

    printf("fport: %02hhx\r\n", get_fport(payload));

    printf("frm payload: ");
    print_frm_payload(payload, payload_len);

    printf("---------\r\n");
    printf("dev addr: %08x (big endian)\r\n", get_dev_addr(payload));
    printf("fctrl: %02hhx\r\n", get_fctrl(payload));
    printf("fcnt: %04x (big endian)\r\n", get_fcnt(payload));

    printf("---------\r\n");
    printf("message type: %s (0x%02hhx)\r\n", get_msg_type_str(get_msg_type(payload)), get_msg_type(payload));
    printf(" message ver: %s (0x%02hhx)\r\n", get_msg_ver_str(get_msg_ver(payload)), get_msg_ver(payload));
    printf("   direction: %s\r\n", get_msg_dir(get_msg_type(payload)));
    printf("        fcnt: %d\r\n", get_fcnt(payload));
}
