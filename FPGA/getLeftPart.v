// ==============================================================
// RTL generated by Vivado(TM) HLS - High-Level Synthesis from C, C++ and OpenCL
// Version: 2020.1
// Copyright (C) 1986-2020 Xilinx, Inc. All Rights Reserved.
// 
// ===========================================================

`timescale 1 ns / 1 ps 

module getLeftPart (
        ap_clk,
        ap_rst,
        ap_start,
        ap_done,
        ap_idle,
        ap_ready,
        slotIndex,
        auxiliaryList_address0,
        auxiliaryList_ce0,
        auxiliaryList_q0,
        ap_return_0,
        ap_return_1,
        ap_return_2
);

parameter    ap_ST_fsm_state1 = 4'd1;
parameter    ap_ST_fsm_state2 = 4'd2;
parameter    ap_ST_fsm_state3 = 4'd4;
parameter    ap_ST_fsm_state4 = 4'd8;

input   ap_clk;
input   ap_rst;
input   ap_start;
output   ap_done;
output   ap_idle;
output   ap_ready;
input  [13:0] slotIndex;
output  [13:0] auxiliaryList_address0;
output   auxiliaryList_ce0;
input  [63:0] auxiliaryList_q0;
output  [7:0] ap_return_0;
output  [63:0] ap_return_1;
output  [63:0] ap_return_2;

reg ap_done;
reg ap_idle;
reg ap_ready;
reg auxiliaryList_ce0;
reg[7:0] ap_return_0;
reg[63:0] ap_return_1;
reg[63:0] ap_return_2;

(* fsm_encoding = "none" *) reg   [3:0] ap_CS_fsm;
wire    ap_CS_fsm_state1;
wire   [19:0] mul_ln45_fu_558_p2;
wire   [5:0] slotBitIdxInWord_fu_194_p1;
reg   [5:0] slotBitIdxInWord_reg_569;
wire    ap_CS_fsm_state2;
wire   [2:0] temp_fu_208_p2;
reg   [2:0] temp_reg_577;
wire   [6:0] select_ln64_fu_268_p3;
reg   [6:0] select_ln64_reg_582;
wire   [0:0] icmp_ln53_fu_202_p2;
wire   [19:0] bitIdx_fu_299_p2;
reg   [19:0] bitIdx_reg_593;
wire   [8:0] extractedBitsNum_fu_313_p2;
reg   [8:0] extractedBitsNum_reg_598;
wire   [5:0] LPBitInWord_fu_319_p2;
reg   [5:0] LPBitInWord_reg_603;
reg   [2:0] trunc_ln8_reg_609;
wire   [0:0] icmp_ln90_fu_335_p2;
reg   [0:0] icmp_ln90_reg_615;
wire   [63:0] select_ln80_1_fu_450_p3;
wire    ap_CS_fsm_state3;
wire   [63:0] select_ln80_2_fu_458_p3;
reg   [63:0] leftPart_0_reg_71;
reg   [63:0] leftPart12_0_reg_83;
reg   [2:0] temp_0_reg_95;
reg   [5:0] LPBitInWord_0_reg_106;
reg   [2:0] LPWordIdx_0_reg_118;
reg   [8:0] extractedBitsNum_0_reg_130;
reg   [19:0] slotBitIdxInWord_0_i_reg_141;
reg   [63:0] leftPart_4_reg_150;
reg   [63:0] leftPart12_4_reg_160;
reg   [5:0] LPBitInWord_1_reg_170;
reg   [2:0] LPWordIdx_1_reg_180;
wire   [63:0] zext_ln71_fu_294_p1;
wire   [6:0] zext_ln47_1_fu_220_p1;
wire   [6:0] tempB_fu_224_p2;
wire   [8:0] toExtractBitsNum_7_fu_214_p2;
wire   [8:0] zext_ln56_fu_230_p1;
wire   [0:0] icmp_ln57_fu_234_p2;
wire   [8:0] toExtractBitsNum_5_fu_240_p3;
wire   [6:0] zext_ln47_fu_198_p1;
wire   [6:0] toExtractBitsNum_fu_252_p2;
wire   [8:0] zext_ln65_fu_258_p1;
wire   [0:0] icmp_ln64_fu_262_p2;
wire   [6:0] trunc_ln57_fu_248_p1;
wire   [13:0] slotWordIdx2_fu_284_p4;
wire   [19:0] zext_ln64_1_fu_280_p1;
wire   [8:0] zext_ln64_fu_276_p1;
wire   [5:0] trunc_ln86_fu_305_p1;
wire   [5:0] trunc_ln86_1_fu_309_p1;
wire   [63:0] zext_ln74_fu_346_p1;
wire   [63:0] shl_ln74_fu_349_p2;
wire   [63:0] zext_ln75_fu_361_p1;
wire   [63:0] lshr_ln75_fu_364_p2;
wire   [63:0] extractPartMask_fu_355_p2;
wire   [0:0] icmp_ln70_fu_341_p2;
wire   [63:0] extractPart_2_fu_370_p2;
wire   [0:0] trunc_ln78_fu_390_p1;
wire   [0:0] icmp_ln77_fu_384_p2;
wire   [63:0] select_ln78_fu_394_p3;
wire   [63:0] select_ln78_1_fu_402_p3;
wire   [63:0] extractPart_3_fu_376_p3;
wire   [63:0] zext_ln80_fu_426_p1;
wire   [63:0] select_ln77_1_fu_418_p3;
wire   [63:0] select_ln77_fu_410_p3;
wire   [63:0] select_ln80_fu_436_p3;
wire   [63:0] shl_ln80_fu_430_p2;
wire   [63:0] add_ln80_fu_444_p2;
wire    ap_CS_fsm_state4;
wire   [0:0] trunc_ln95_fu_470_p1;
wire   [6:0] zext_ln95_1_fu_466_p1;
wire   [6:0] add_ln95_fu_482_p2;
wire  signed [31:0] sext_ln95_fu_488_p1;
wire   [63:0] select_ln95_fu_474_p3;
wire   [63:0] zext_ln95_fu_492_p1;
wire   [63:0] lshr_ln95_fu_496_p2;
wire   [63:0] shl_ln96_fu_506_p2;
wire   [63:0] xor_ln96_fu_512_p2;
wire   [63:0] and_ln96_fu_518_p2;
wire   [7:0] counter_fu_502_p1;
wire   [63:0] select_ln96_fu_524_p3;
wire   [63:0] select_ln96_1_fu_532_p3;
wire   [13:0] mul_ln45_fu_558_p0;
wire   [7:0] mul_ln45_fu_558_p1;
reg   [7:0] ap_return_0_preg;
reg   [63:0] ap_return_1_preg;
reg   [63:0] ap_return_2_preg;
reg   [3:0] ap_NS_fsm;
wire   [19:0] mul_ln45_fu_558_p00;

// power-on initialization
initial begin
#0 ap_CS_fsm = 4'd1;
#0 ap_return_0_preg = 8'd0;
#0 ap_return_1_preg = 64'd0;
#0 ap_return_2_preg = 64'd0;
end

insert_mul_mul_14bkb #(
    .ID( 1 ),
    .NUM_STAGE( 1 ),
    .din0_WIDTH( 14 ),
    .din1_WIDTH( 8 ),
    .dout_WIDTH( 20 ))
insert_mul_mul_14bkb_U4(
    .din0(mul_ln45_fu_558_p0),
    .din1(mul_ln45_fu_558_p1),
    .dout(mul_ln45_fu_558_p2)
);

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_CS_fsm <= ap_ST_fsm_state1;
    end else begin
        ap_CS_fsm <= ap_NS_fsm;
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_return_0_preg <= 8'd0;
    end else begin
        if ((1'b1 == ap_CS_fsm_state4)) begin
            ap_return_0_preg <= counter_fu_502_p1;
        end
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_return_1_preg <= 64'd0;
    end else begin
        if ((1'b1 == ap_CS_fsm_state4)) begin
            ap_return_1_preg <= select_ln96_fu_524_p3;
        end
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_return_2_preg <= 64'd0;
    end else begin
        if ((1'b1 == ap_CS_fsm_state4)) begin
            ap_return_2_preg <= select_ln96_1_fu_532_p3;
        end
    end
end

always @ (posedge ap_clk) begin
    if (((icmp_ln90_reg_615 == 1'd0) & (1'b1 == ap_CS_fsm_state3))) begin
        LPBitInWord_0_reg_106 <= LPBitInWord_reg_603;
    end else if (((ap_start == 1'b1) & (1'b1 == ap_CS_fsm_state1))) begin
        LPBitInWord_0_reg_106 <= 6'd0;
    end
end

always @ (posedge ap_clk) begin
    if (((icmp_ln90_reg_615 == 1'd1) & (1'b1 == ap_CS_fsm_state3))) begin
        LPBitInWord_1_reg_170 <= LPBitInWord_reg_603;
    end else if (((icmp_ln53_fu_202_p2 == 1'd1) & (1'b1 == ap_CS_fsm_state2))) begin
        LPBitInWord_1_reg_170 <= LPBitInWord_0_reg_106;
    end
end

always @ (posedge ap_clk) begin
    if (((icmp_ln90_reg_615 == 1'd0) & (1'b1 == ap_CS_fsm_state3))) begin
        LPWordIdx_0_reg_118 <= trunc_ln8_reg_609;
    end else if (((ap_start == 1'b1) & (1'b1 == ap_CS_fsm_state1))) begin
        LPWordIdx_0_reg_118 <= 3'd0;
    end
end

always @ (posedge ap_clk) begin
    if (((icmp_ln90_reg_615 == 1'd1) & (1'b1 == ap_CS_fsm_state3))) begin
        LPWordIdx_1_reg_180 <= trunc_ln8_reg_609;
    end else if (((icmp_ln53_fu_202_p2 == 1'd1) & (1'b1 == ap_CS_fsm_state2))) begin
        LPWordIdx_1_reg_180 <= LPWordIdx_0_reg_118;
    end
end

always @ (posedge ap_clk) begin
    if (((icmp_ln90_reg_615 == 1'd0) & (1'b1 == ap_CS_fsm_state3))) begin
        extractedBitsNum_0_reg_130 <= extractedBitsNum_reg_598;
    end else if (((ap_start == 1'b1) & (1'b1 == ap_CS_fsm_state1))) begin
        extractedBitsNum_0_reg_130 <= 9'd0;
    end
end

always @ (posedge ap_clk) begin
    if (((icmp_ln90_reg_615 == 1'd0) & (1'b1 == ap_CS_fsm_state3))) begin
        leftPart12_0_reg_83 <= select_ln80_2_fu_458_p3;
    end else if (((ap_start == 1'b1) & (1'b1 == ap_CS_fsm_state1))) begin
        leftPart12_0_reg_83 <= 64'd0;
    end
end

always @ (posedge ap_clk) begin
    if (((icmp_ln90_reg_615 == 1'd1) & (1'b1 == ap_CS_fsm_state3))) begin
        leftPart12_4_reg_160 <= select_ln80_2_fu_458_p3;
    end else if (((icmp_ln53_fu_202_p2 == 1'd1) & (1'b1 == ap_CS_fsm_state2))) begin
        leftPart12_4_reg_160 <= leftPart12_0_reg_83;
    end
end

always @ (posedge ap_clk) begin
    if (((icmp_ln90_reg_615 == 1'd0) & (1'b1 == ap_CS_fsm_state3))) begin
        leftPart_0_reg_71 <= select_ln80_1_fu_450_p3;
    end else if (((ap_start == 1'b1) & (1'b1 == ap_CS_fsm_state1))) begin
        leftPart_0_reg_71 <= 64'd0;
    end
end

always @ (posedge ap_clk) begin
    if (((icmp_ln90_reg_615 == 1'd1) & (1'b1 == ap_CS_fsm_state3))) begin
        leftPart_4_reg_150 <= select_ln80_1_fu_450_p3;
    end else if (((icmp_ln53_fu_202_p2 == 1'd1) & (1'b1 == ap_CS_fsm_state2))) begin
        leftPart_4_reg_150 <= leftPart_0_reg_71;
    end
end

always @ (posedge ap_clk) begin
    if (((icmp_ln90_reg_615 == 1'd0) & (1'b1 == ap_CS_fsm_state3))) begin
        slotBitIdxInWord_0_i_reg_141 <= bitIdx_reg_593;
    end else if (((ap_start == 1'b1) & (1'b1 == ap_CS_fsm_state1))) begin
        slotBitIdxInWord_0_i_reg_141 <= mul_ln45_fu_558_p2;
    end
end

always @ (posedge ap_clk) begin
    if (((icmp_ln90_reg_615 == 1'd0) & (1'b1 == ap_CS_fsm_state3))) begin
        temp_0_reg_95 <= temp_reg_577;
    end else if (((ap_start == 1'b1) & (1'b1 == ap_CS_fsm_state1))) begin
        temp_0_reg_95 <= 3'd0;
    end
end

always @ (posedge ap_clk) begin
    if (((icmp_ln53_fu_202_p2 == 1'd0) & (1'b1 == ap_CS_fsm_state2))) begin
        LPBitInWord_reg_603 <= LPBitInWord_fu_319_p2;
        bitIdx_reg_593 <= bitIdx_fu_299_p2;
        extractedBitsNum_reg_598 <= extractedBitsNum_fu_313_p2;
        icmp_ln90_reg_615 <= icmp_ln90_fu_335_p2;
        select_ln64_reg_582 <= select_ln64_fu_268_p3;
        trunc_ln8_reg_609 <= {{extractedBitsNum_fu_313_p2[8:6]}};
    end
end

always @ (posedge ap_clk) begin
    if ((1'b1 == ap_CS_fsm_state2)) begin
        slotBitIdxInWord_reg_569 <= slotBitIdxInWord_fu_194_p1;
        temp_reg_577 <= temp_fu_208_p2;
    end
end

always @ (*) begin
    if (((1'b1 == ap_CS_fsm_state4) | ((ap_start == 1'b0) & (1'b1 == ap_CS_fsm_state1)))) begin
        ap_done = 1'b1;
    end else begin
        ap_done = 1'b0;
    end
end

always @ (*) begin
    if (((ap_start == 1'b0) & (1'b1 == ap_CS_fsm_state1))) begin
        ap_idle = 1'b1;
    end else begin
        ap_idle = 1'b0;
    end
end

always @ (*) begin
    if ((1'b1 == ap_CS_fsm_state4)) begin
        ap_ready = 1'b1;
    end else begin
        ap_ready = 1'b0;
    end
end

always @ (*) begin
    if ((1'b1 == ap_CS_fsm_state4)) begin
        ap_return_0 = counter_fu_502_p1;
    end else begin
        ap_return_0 = ap_return_0_preg;
    end
end

always @ (*) begin
    if ((1'b1 == ap_CS_fsm_state4)) begin
        ap_return_1 = select_ln96_fu_524_p3;
    end else begin
        ap_return_1 = ap_return_1_preg;
    end
end

always @ (*) begin
    if ((1'b1 == ap_CS_fsm_state4)) begin
        ap_return_2 = select_ln96_1_fu_532_p3;
    end else begin
        ap_return_2 = ap_return_2_preg;
    end
end

always @ (*) begin
    if ((1'b1 == ap_CS_fsm_state2)) begin
        auxiliaryList_ce0 = 1'b1;
    end else begin
        auxiliaryList_ce0 = 1'b0;
    end
end

always @ (*) begin
    case (ap_CS_fsm)
        ap_ST_fsm_state1 : begin
            if (((ap_start == 1'b1) & (1'b1 == ap_CS_fsm_state1))) begin
                ap_NS_fsm = ap_ST_fsm_state2;
            end else begin
                ap_NS_fsm = ap_ST_fsm_state1;
            end
        end
        ap_ST_fsm_state2 : begin
            if (((icmp_ln53_fu_202_p2 == 1'd1) & (1'b1 == ap_CS_fsm_state2))) begin
                ap_NS_fsm = ap_ST_fsm_state4;
            end else begin
                ap_NS_fsm = ap_ST_fsm_state3;
            end
        end
        ap_ST_fsm_state3 : begin
            if (((icmp_ln90_reg_615 == 1'd1) & (1'b1 == ap_CS_fsm_state3))) begin
                ap_NS_fsm = ap_ST_fsm_state4;
            end else begin
                ap_NS_fsm = ap_ST_fsm_state2;
            end
        end
        ap_ST_fsm_state4 : begin
            ap_NS_fsm = ap_ST_fsm_state1;
        end
        default : begin
            ap_NS_fsm = 'bx;
        end
    endcase
end

assign LPBitInWord_fu_319_p2 = (trunc_ln86_fu_305_p1 + trunc_ln86_1_fu_309_p1);

assign add_ln80_fu_444_p2 = (select_ln80_fu_436_p3 + shl_ln80_fu_430_p2);

assign add_ln95_fu_482_p2 = ($signed(7'd126) + $signed(zext_ln95_1_fu_466_p1));

assign and_ln96_fu_518_p2 = (xor_ln96_fu_512_p2 & select_ln95_fu_474_p3);

assign ap_CS_fsm_state1 = ap_CS_fsm[32'd0];

assign ap_CS_fsm_state2 = ap_CS_fsm[32'd1];

assign ap_CS_fsm_state3 = ap_CS_fsm[32'd2];

assign ap_CS_fsm_state4 = ap_CS_fsm[32'd3];

assign auxiliaryList_address0 = zext_ln71_fu_294_p1;

assign bitIdx_fu_299_p2 = (zext_ln64_1_fu_280_p1 + slotBitIdxInWord_0_i_reg_141);

assign counter_fu_502_p1 = lshr_ln95_fu_496_p2[7:0];

assign extractPartMask_fu_355_p2 = ($signed(64'd18446744073709551615) + $signed(shl_ln74_fu_349_p2));

assign extractPart_2_fu_370_p2 = (lshr_ln75_fu_364_p2 & extractPartMask_fu_355_p2);

assign extractPart_3_fu_376_p3 = ((icmp_ln70_fu_341_p2[0:0] === 1'b1) ? auxiliaryList_q0 : extractPart_2_fu_370_p2);

assign extractedBitsNum_fu_313_p2 = (zext_ln64_fu_276_p1 + extractedBitsNum_0_reg_130);

assign icmp_ln53_fu_202_p2 = ((temp_0_reg_95 == 3'd4) ? 1'b1 : 1'b0);

assign icmp_ln57_fu_234_p2 = ((toExtractBitsNum_7_fu_214_p2 < zext_ln56_fu_230_p1) ? 1'b1 : 1'b0);

assign icmp_ln64_fu_262_p2 = ((toExtractBitsNum_5_fu_240_p3 > zext_ln65_fu_258_p1) ? 1'b1 : 1'b0);

assign icmp_ln70_fu_341_p2 = ((select_ln64_reg_582 == 7'd64) ? 1'b1 : 1'b0);

assign icmp_ln77_fu_384_p2 = ((LPBitInWord_0_reg_106 == 6'd0) ? 1'b1 : 1'b0);

assign icmp_ln90_fu_335_p2 = ((extractedBitsNum_fu_313_p2 > 9'd78) ? 1'b1 : 1'b0);

assign lshr_ln75_fu_364_p2 = auxiliaryList_q0 >> zext_ln75_fu_361_p1;

assign lshr_ln95_fu_496_p2 = select_ln95_fu_474_p3 >> zext_ln95_fu_492_p1;

assign mul_ln45_fu_558_p0 = mul_ln45_fu_558_p00;

assign mul_ln45_fu_558_p00 = slotIndex;

assign mul_ln45_fu_558_p1 = 20'd79;

assign select_ln64_fu_268_p3 = ((icmp_ln64_fu_262_p2[0:0] === 1'b1) ? toExtractBitsNum_fu_252_p2 : trunc_ln57_fu_248_p1);

assign select_ln77_1_fu_418_p3 = ((icmp_ln77_fu_384_p2[0:0] === 1'b1) ? select_ln78_1_fu_402_p3 : leftPart12_0_reg_83);

assign select_ln77_fu_410_p3 = ((icmp_ln77_fu_384_p2[0:0] === 1'b1) ? select_ln78_fu_394_p3 : leftPart_0_reg_71);

assign select_ln78_1_fu_402_p3 = ((trunc_ln78_fu_390_p1[0:0] === 1'b1) ? 64'd0 : leftPart12_0_reg_83);

assign select_ln78_fu_394_p3 = ((trunc_ln78_fu_390_p1[0:0] === 1'b1) ? leftPart_0_reg_71 : 64'd0);

assign select_ln80_1_fu_450_p3 = ((trunc_ln78_fu_390_p1[0:0] === 1'b1) ? select_ln77_fu_410_p3 : add_ln80_fu_444_p2);

assign select_ln80_2_fu_458_p3 = ((trunc_ln78_fu_390_p1[0:0] === 1'b1) ? add_ln80_fu_444_p2 : select_ln77_1_fu_418_p3);

assign select_ln80_fu_436_p3 = ((trunc_ln78_fu_390_p1[0:0] === 1'b1) ? select_ln77_1_fu_418_p3 : select_ln77_fu_410_p3);

assign select_ln95_fu_474_p3 = ((trunc_ln95_fu_470_p1[0:0] === 1'b1) ? leftPart12_4_reg_160 : leftPart_4_reg_150);

assign select_ln96_1_fu_532_p3 = ((trunc_ln95_fu_470_p1[0:0] === 1'b1) ? and_ln96_fu_518_p2 : leftPart12_4_reg_160);

assign select_ln96_fu_524_p3 = ((trunc_ln95_fu_470_p1[0:0] === 1'b1) ? leftPart_4_reg_150 : and_ln96_fu_518_p2);

assign sext_ln95_fu_488_p1 = $signed(add_ln95_fu_482_p2);

assign shl_ln74_fu_349_p2 = 64'd1 << zext_ln74_fu_346_p1;

assign shl_ln80_fu_430_p2 = extractPart_3_fu_376_p3 << zext_ln80_fu_426_p1;

assign shl_ln96_fu_506_p2 = 64'd3 << zext_ln95_fu_492_p1;

assign slotBitIdxInWord_fu_194_p1 = slotBitIdxInWord_0_i_reg_141[5:0];

assign slotWordIdx2_fu_284_p4 = {{slotBitIdxInWord_0_i_reg_141[19:6]}};

assign tempB_fu_224_p2 = ($signed(7'd64) - $signed(zext_ln47_1_fu_220_p1));

assign temp_fu_208_p2 = (3'd1 + temp_0_reg_95);

assign toExtractBitsNum_5_fu_240_p3 = ((icmp_ln57_fu_234_p2[0:0] === 1'b1) ? toExtractBitsNum_7_fu_214_p2 : zext_ln56_fu_230_p1);

assign toExtractBitsNum_7_fu_214_p2 = (9'd79 - extractedBitsNum_0_reg_130);

assign toExtractBitsNum_fu_252_p2 = ($signed(7'd64) - $signed(zext_ln47_fu_198_p1));

assign trunc_ln57_fu_248_p1 = toExtractBitsNum_5_fu_240_p3[6:0];

assign trunc_ln78_fu_390_p1 = LPWordIdx_0_reg_118[0:0];

assign trunc_ln86_1_fu_309_p1 = select_ln64_fu_268_p3[5:0];

assign trunc_ln86_fu_305_p1 = extractedBitsNum_0_reg_130[5:0];

assign trunc_ln95_fu_470_p1 = LPWordIdx_1_reg_180[0:0];

assign xor_ln96_fu_512_p2 = (shl_ln96_fu_506_p2 ^ 64'd18446744073709551615);

assign zext_ln47_1_fu_220_p1 = LPBitInWord_0_reg_106;

assign zext_ln47_fu_198_p1 = slotBitIdxInWord_fu_194_p1;

assign zext_ln56_fu_230_p1 = tempB_fu_224_p2;

assign zext_ln64_1_fu_280_p1 = select_ln64_fu_268_p3;

assign zext_ln64_fu_276_p1 = select_ln64_fu_268_p3;

assign zext_ln65_fu_258_p1 = toExtractBitsNum_fu_252_p2;

assign zext_ln71_fu_294_p1 = slotWordIdx2_fu_284_p4;

assign zext_ln74_fu_346_p1 = select_ln64_reg_582;

assign zext_ln75_fu_361_p1 = slotBitIdxInWord_reg_569;

assign zext_ln80_fu_426_p1 = LPBitInWord_0_reg_106;

assign zext_ln95_1_fu_466_p1 = LPBitInWord_1_reg_170;

assign zext_ln95_fu_492_p1 = $unsigned(sext_ln95_fu_488_p1);

endmodule //getLeftPart
