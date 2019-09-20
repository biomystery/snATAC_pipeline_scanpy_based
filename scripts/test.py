
bf = pysam.AlignmentFile(input_bam, 'rb')

with open('cb.txt','w') as f:
    for read in bf:
        cb_name=read.get_tag('CR')
        #barcode=read.get_tag('CR').split('-1')[0]
        barcode=cb_name
        f.write(read.query_name+','+cb_name+','+barcode+'\n')

bf.close()

