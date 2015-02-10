from Uniprot_Parser import *
from seq_class import *
import urllib
record=SeqIO.read("sequence.gb", "genbank")
locus=['NGO0243', 'NGO0244', 'NGO0245', 'NGO_t05', 'NGO0248', 'NGO0249', 'NGO0250', 'NGO0250a', 'NGO0252', 'NGO0253', 'NGO0255', 'NGO0256', 'NGO0257', 'NGO0258', 'NGO0259', 'NGO0260', 'NGO0261', 'NGO0262', 'NGO0263', 'NGO0264', 'NGO0265', 'NGO0266', 'NGO0267', 'NGO0268', 'NGO0269', 'NGO0270', 'NGO0271', 'NGO0272', 'NGO272a', 'NGO0274', 'NGO0275', 'NGO0276', 'NGO0277', 'NGO0278', 'NGO0281', 'NGO0282', 'NGO0283', 'NGO0284', 'NGO0285', 'NGO0286', 'NGO0288', 'NGO0289', 'NGO0290', 'NGO0291', 'NGO0292', 'NGO0293', 'NGO0294', 'NGO_t06', 'NGO0295', 'NGO0296', 'NGO0297', 'NGO0298', 'NGO0299', 'NGO0300', 'NGO0302', 'NGO0303', 'NGO0304', 'NGO0305', 'NGO0306', 'NGO0307', 'NGO0308', 'NGO0309', 'NGO0310', 'NGO0311', 'NGO0312', 'NGO0313', 'NGO0314', 'NGO0315', 'NGO0316', 'NGO0317', 'NGO0318', 'NGO_t07', 'NGO_t08', 'NGO0319', 'NGO0320', 'NGO0321', 'NGO0322', 'NGO0323', 'NGO0324', 'NGO0326', 'NGO0327', 'NGO0328', 'NGO0329', 'NGO0329a', 'NGO0331', 'NGO0332', 'NGO0333', 'NGO0335', 'NGO0336', 'NGO0337', 'NGO0338', 'NGO0339', 'NGO0340', 'NGO0342', 'NGO0343', 'NGO0344', 'NGO0345', 'NGO0346', 'NGO0347', 'NGO0348', 'NGO0349', 'NGO0350', 'NGO0351', 'NGO0352', 'NGO0353', 'NGO0354', 'NGO0355', 'NGO0356', 'NGO0357', 'NGO0358', 'NGO0359', 'NGO0360', 'NGO0361', 'NGO0362', 'NGO0363', 'NGO0364', 'NGO0365', 'NGO0366', 'NGO0367', 'NGO0368', 'NGO0369', 'NGO0370', 'NGO0371', 'NGO0372', 'NGO0373', 'NGO0374', 'NGO0375', 'NGO0376', 'NGO0377', 'NGO0378', 'NGO0379', 'NGO0380', 'NGO0381', 'NGO0382', 'NGO0383', 'NGO0384', 'NGO0385', 'NGO0386', 'NGO0387', 'NGO0388', 'NGO0389', 'NGO0390', 'NGO0391', 'NGO0392', 'NGO0393', 'NGO0394', 'NGO0395', 'NGO0397', 'NGO0398', 'NGO0399', 'NGO0400', 'NGO0401', 'NGO0402', 'NGO0403', 'NGO0404', 'NGO0405', 'NGO0407', 'NGO0408', 'NGO0409', 'NGO0410', 'NGO0411', 'NGO0412', 'NGO0413', 'NGO0414', 'NGO0415', 'NGO0416', 'NGO0417', 'NGO0418', 'NGO0419', 'NGO0420', 'NGO0421', 'NGO0422', 'NGO0423', 'NGO0424', 'NGO0425', 'NGO0426', 'NGO0427', 'NGO0428', 'NGO0429', 'NGO0430', 'NGO0432', 'NGO0433', 'NGO0434', 'NGO0435', 'NGO_t09', 'NGO0436', 'NGO0437', 'NGO0439', 'NGO0440', 'NGO_t10', 'NGO_t11', 'NGO0441', 'NGO0442', 'NGO0443', 'NGO0444', 'NGO0445', 'NGO0446', 'NGO0448', 'NGO0449', 'NGO0451', 'NGO0452', 'NGO0453', 'NGO0454', 'NGO0455', 'NGO0456', 'NGO0457', 'NGO0459', 'NGO0460', 'NGO0461', 'NGO_t12', 'NGO_t13', 'NGO0462', 'NGO0463', 'NGO0464', 'NGO0465', 'NGO0466', 'NGO0467', 'NGO0468', 'NGO0469', 'NGO0470', 'NGO0472', 'NGO0473', 'NGO0474', 'NGO0475', 'NGO0476', 'NGO0478', 'NGO0479', 'NGO0480', 'NGO0481', 'NGO0483', 'NGO0484', 'NGO0485']
proteins=['YP_207408.1', 'YP_207409.1', 'YP_207410.1', 'YP_207413.1', 'YP_207414.1', 'YP_207415.1', 'YP_008914848.1', 'YP_207416.1', 'YP_207417.1', 'YP_207418.1', 'YP_207419.1', 'YP_207420.1', 'YP_207421.1', 'YP_207422.1', 'YP_207423.1', 'YP_207424.1', 'YP_207425.1', 'YP_207426.1', 'YP_207427.1', 'YP_207428.1', 'YP_207429.1', 'YP_207430.1', 'YP_207431.1', 'YP_207432.1', 'YP_207433.1', 'YP_207434.1', 'YP_207435.1', 'YP_009115479.1', 'YP_207436.1', 'YP_207437.1', 'YP_207438.1', 'YP_207439.1', 'YP_207440.1', 'YP_207442.1', 'YP_207443.1', 'YP_207444.1', 'YP_207445.1', 'YP_207446.1', 'YP_207447.1', 'YP_207448.1', 'YP_207449.1', 'YP_207451.1', 'YP_207452.1', 'YP_207453.1', 'YP_207454.1', 'YP_207455.1', 'YP_207456.1', 'YP_207457.1', 'YP_207458.1', 'YP_207459.2', 'YP_207460.1', 'YP_207461.1', 'YP_207462.1', 'YP_207463.1', 'YP_207464.1', 'YP_207465.1', 'YP_207466.1', 'YP_207467.1', 'YP_207468.1', 'YP_207469.1', 'YP_207470.1', 'YP_207471.1', 'YP_207472.1', 'YP_207473.1', 'YP_207474.1', 'YP_207475.1', 'YP_207476.1', 'YP_207477.1', 'YP_207478.1', 'YP_207479.1', 'YP_207480.1', 'YP_207481.1', 'YP_207482.1', 'YP_207483.1', 'YP_207484.1', 'YP_207485.1', 'YP_207486.1', 'YP_207487.1', 'YP_009115480.1', 'YP_207489.1', 'YP_207490.1', 'YP_207491.1', 'YP_207492.1', 'YP_207493.1', 'YP_207494.1', 'YP_207495.1', 'YP_207496.1', 'YP_207497.1', 'YP_207498.1', 'YP_207499.1', 'YP_207500.1', 'YP_207501.1', 'YP_207502.1', 'YP_207503.1', 'YP_207504.1', 'YP_207505.1', 'YP_207506.1', 'YP_207507.1', 'YP_207508.1', 'YP_207509.1', 'YP_207510.1', 'YP_207511.1', 'YP_207512.1', 'YP_207513.1', 'YP_207514.1', 'YP_207515.1', 'YP_207516.1', 'YP_207517.1', 'YP_207518.1', 'YP_207519.1', 'YP_207520.1', 'YP_207521.1', 'YP_207522.1', 'YP_207523.1', 'YP_207524.2', 'YP_207525.2', 'YP_207526.1', 'YP_207527.1', 'YP_207528.1', 'YP_207529.1', 'YP_207530.1', 'YP_207531.1', 'YP_207532.1', 'YP_207533.1', 'YP_207534.1', 'YP_207535.1', 'YP_207536.1', 'YP_207537.1', 'YP_207538.1', 'YP_207539.1', 'YP_207540.1', 'YP_207541.1', 'YP_207542.1', 'YP_207543.1', 'YP_207544.1', 'YP_207545.1', 'YP_207546.1', 'YP_207547.1', 'YP_207548.1', 'YP_207549.1', 'YP_207550.1', 'YP_207551.1', 'YP_207553.1', 'YP_207554.1', 'YP_207555.1', 'YP_207556.1', 'YP_207557.1', 'YP_207558.1', 'YP_207559.1', 'YP_207560.1', 'YP_207561.1', 'YP_207562.1', 'YP_207563.1', 'YP_207564.1', 'YP_207565.1', 'YP_207566.1', 'YP_207567.1', 'YP_207568.1', 'YP_207569.1', 'YP_207570.1', 'YP_207571.1', 'YP_207572.1', 'YP_207573.1', 'YP_207574.1', 'YP_207575.1', 'YP_207576.1', 'YP_207578.1', 'YP_207579.1', 'YP_207580.1', 'YP_207581.1', 'YP_207582.1', 'YP_207583.1', 'YP_207584.1', 'YP_207585.1', 'YP_207586.1', 'YP_207587.1', 'YP_207588.1', 'YP_207589.1', 'YP_207590.1', 'YP_207591.1', 'YP_207592.1', 'YP_207593.1', 'YP_207594.1', 'YP_207595.1', 'YP_207596.1', 'YP_207597.1', 'YP_207598.1', 'YP_207599.1', 'YP_207600.1', 'YP_207601.1', 'YP_207602.1', 'YP_207603.1', 'YP_207604.1', 'YP_207605.1', 'YP_207606.1', 'YP_207607.1', 'YP_207608.1', 'YP_207609.1', 'YP_207610.1', 'YP_207611.1', 'YP_207612.1', 'YP_207613.1', 'YP_207614.1', 'YP_207615.1', 'YP_207616.1', 'YP_207617.1', 'YP_207618.1', 'YP_207619.1', 'YP_207620.1', 'YP_207622.1', 'YP_207623.1', 'YP_207624.1', 'YP_207625.1', 'YP_207626.1', 'YP_207628.1', 'YP_207629.1', 'YP_207630.1', 'YP_207631.1', 'YP_207633.1', 'YP_207634.1', 'YP_207635.1']




def uniprot_ID(record,proteins):
    IDs=[]
    for j in range(len(proteins)):
        data = urllib.request.urlopen("http://www.uniprot.org/uniprot/?query="+proteins[j]+"&sort=score").read()
        if len(data.split())>=3312:
            x=str(data.split()[3312])
            IDs.append(x[6:12])
        else:
            IDs.append(proteins[j])
    return IDs



def info_uniprot():
    identifier=[]
    handle = open("../res/"+"uniprot.txt")
    records = parse(handle) # Uses the function 'parse' from the module. 
    for record in records:
        for i in range(len(ID)):
            identifier.append([])
            if record["AC"]==ID[i]+";":
                #identifier.append([])
                identifier[i].append(ID[i])
                identifier[i].append(record["ID"])
                identifier[i].append(record["DE"])
                identifier[i].append(record["CC"])
    handle.close()
    return identifier


def more_info_uniprot():
    handle = open("../res/"+"uniprot_XML.xml")
    records=UniprotIO.UniprotIterator(handle,return_raw_comments=True)
    refs=[]
    for record in records:
        for i in range(len(ID)):
            if record.id==ID[i]:
                refs.append([ID[i]]+record.dbxrefs)#GO´s
    handle.close()
    return refs
    
print(uniprot_ID(record,proteins))

#data = urllib.request.urlopen("http://www.uniprot.org/uniprot/?query="+"YP_207409.1"+"&sort=score").read()
#print(data.split().index(b'id="Q5F9Z0"'))