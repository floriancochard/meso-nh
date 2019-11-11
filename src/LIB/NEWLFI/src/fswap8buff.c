/*
*MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
*MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
*MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
*MNH_LIC for details. version 1.
*/#define SWAP4(w) ((w<<24)|((w<<8)&0xff0000)|((w>>8)&0xff00)|((w>>24)&0xff))

void swap8_(int *outbuff,int *inbuff,int *nbint8)
{
    int i;

    for (i=0; i<*nbint8; i++){
	outbuff[0] = SWAP4(inbuff[1]);
	outbuff[1] = SWAP4(inbuff[0]);
	outbuff += 2;
	inbuff  += 2;
    }
}
