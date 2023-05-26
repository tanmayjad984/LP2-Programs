LP2 


IS

#####1) ) Write Java/Python/C++ code to accept a string from the User. Perform bit-wise AND and
XOR each character in the string with 127 and display the results. #####


#Subs_Cipher.java

public class Subs_Cipher {
public static void main(String[] args) {
String str = "HelloWorld";
StringBuilder andResult = new StringBuilder();
StringBuilder xorResult = new StringBuilder();
for (char c : str.toCharArray()) {
// Perform AND operation
char andChar = (char) (c & 127);
andResult.append(andChar);
// Perform XOR operation
char xorChar = (char) (c ^ 127);
xorResult.append(xorChar);
}
// Display results
System.out.println("XOR Result: " + xorResult.toString());
System.out.println("Original String: " + str);
System.out.println("AND Result: " + andResult.toString());
}
}



##### 2) Write Java/Python/C++ code to accept a string from the User. Perform Encryption and
Decryption of the text using any one Transposition technique.  ####

#RowTranspositionCipher.java



import java.util.*;
public class RowTranspositionCipher {
public static void main(String[] args) {
System.out.println("Enter the Plaintext");
Scanner sc = new Scanner(System.in);
String plaintext = sc.next();
System.out.println("Enter the Key");
String key = sc.next();
String ciphertext = encrypt(plaintext, key);
System.out.println("Ciphertext: " + ciphertext);
String decryptedText = decrypt(ciphertext, key);
System.out.println("Decrypted text: " + decryptedText);
}
public static String encrypt(String plaintext, String key) {
// Remove spaces and convert to upper case
plaintext = plaintext.replaceAll("\\s+", "").toUpperCase();
int keyLength = key.length();
int textLength = plaintext.length();
// Calculate number of rows needed to fit the plaintext
int numRows = (int) Math.ceil((double)textLength / keyLength);
// Create a 2D character array to store the plaintext
char[][] matrix = new char[numRows][keyLength];
int k = 0;
for (int i = 0; i < numRows; i++) {
for (int j = 0; j < keyLength; j++) {
if (k < textLength) {
matrix[i][j] = plaintext.charAt(k++);
} else {
matrix[i][j] = 'X'; // Add padding character
}
}
}
// Create an array to hold the ciphertext
char[] ciphertext = new char[textLength];
k = 0;
for (int i = 0; i < keyLength; i++) {
int index = Character.getNumericValue(key.charAt(i)) - 1;
for (int j = 0; j < numRows; j++) {
ciphertext[k++] = matrix[j][index];
}
}
return new String(ciphertext);
}
public static String decrypt(String ciphertext, String key) {
int keyLength = key.length();
int textLength = ciphertext.length();
// Calculate number of rows needed to fit the ciphertext
int numRows = (int) Math.ceil((double)textLength / keyLength);
// Create a 2D character array to store the ciphertext
char[][] matrix = new char[numRows][keyLength];
int k = 0;
for (int i = 0; i < keyLength; i++) {
int index = Character.getNumericValue(key.charAt(i)) - 1;
for (int j = 0; j < numRows; j++) {
matrix[j][index] = ciphertext.charAt(k++);
}
}
// Create an array to hold the plaintext
char[] plaintext = new char[textLength];
k = 0;
for (int i = 0; i < numRows; i++) {
for (int j = 0; j < keyLength; j++) {
if (matrix[i][j] != 'X') {
plaintext[k++] = matrix[i][j];
}
}
}
return new String(plaintext);
}
}



##### 3) Write Java/Python/C++ code to implement DES Encryption Algorithm.   #####



#DesEncryption.java






import java.util.*;
import javax.crypto.Cipher;
import javax.crypto.KeyGenerator;
import javax.crypto.SecretKey;
import javax.crypto.spec.SecretKeySpec;
import java.security.SecureRandom;
public class DesEncryption {
public static void main(String[] args) throws Exception {
// Generate a secret key for DES algorithm
KeyGenerator keyGenerator = KeyGenerator.getInstance("DES");
SecureRandom secureRandom = new SecureRandom();
keyGenerator.init(secureRandom);
SecretKey secretKey = keyGenerator.generateKey();
// Convert the secret key to bytes
byte[] keyBytes = secretKey.getEncoded();
// Create a SecretKeySpec object from the key bytes
SecretKeySpec secretKeySpec = new SecretKeySpec(keyBytes, "DES");
// Create a Cipher object and initialize it with the secret key
Cipher cipher = Cipher.getInstance("DES/ECB/PKCS5Padding");
cipher.init(Cipher.ENCRYPT_MODE, secretKeySpec);
// Encrypt the message
Scanner sc = new Scanner(System.in);
System.out.println("Enter the plaintext");
String message = sc.nextLine();
byte[] encryptedMessageBytes = cipher.doFinal(message.getBytes());
// Print the encrypted message
System.out.println("Encrypted Message: " + new String(encryptedMessageBytes));
// Initialize the Cipher object for decryption
cipher.init(Cipher.DECRYPT_MODE, secretKeySpec);
// Decrypt the message
byte[] decryptedMessageBytes = cipher.doFinal(encryptedMessageBytes);
// Print the decrypted message
System.out.println("Decrypted Message: " + new String(decryptedMessageBytes));
}
}





##### 4) Write Java/Python/C++ code to implement AES Encryption Algorithm.   #####


#AESExample.java


import javax.crypto.Cipher;
import javax.crypto.spec.SecretKeySpec;
import java.nio.charset.StandardCharsets;
import java.util.Base64;
public class AESExample {
private static final String AES_ALGORITHM = "AES";
private static final int KEY_LENGTH = 128;
public static String encrypt(String plaintext, String secretKey) throws Exception {
byte[] keyBytes = fixKeyLength(secretKey.getBytes(StandardCharsets.UTF_8));
SecretKeySpec keySpec = new SecretKeySpec(keyBytes, AES_ALGORITHM);
Cipher cipher = Cipher.getInstance(AES_ALGORITHM);
cipher.init(Cipher.ENCRYPT_MODE, keySpec);
byte[] encryptedBytes = cipher.doFinal(plaintext.getBytes(StandardCharsets.UTF_8));
return Base64.getEncoder().encodeToString(encryptedBytes);
}
public static String decrypt(String ciphertext, String secretKey) throws Exception {
byte[] keyBytes = fixKeyLength(secretKey.getBytes(StandardCharsets.UTF_8));
SecretKeySpec keySpec = new SecretKeySpec(keyBytes, AES_ALGORITHM);
Cipher cipher = Cipher.getInstance(AES_ALGORITHM);
cipher.init(Cipher.DECRYPT_MODE, keySpec);
byte[] decodedBytes = Base64.getDecoder().decode(ciphertext);
byte[] decryptedBytes = cipher.doFinal(decodedBytes);
return new String(decryptedBytes, StandardCharsets.UTF_8);
}
private static byte[] fixKeyLength(byte[] keyBytes) {
int validLength = (KEY_LENGTH / 8);
if (keyBytes.length == validLength) {
return keyBytes;
} else if (keyBytes.length < validLength) {
byte[] newKeyBytes = new byte[validLength];
System.arraycopy(keyBytes, 0, newKeyBytes, 0, keyBytes.length);
return newKeyBytes;
} else {
byte[] newKeyBytes = new byte[validLength];
System.arraycopy(keyBytes, 0, newKeyBytes, 0, validLength);
return newKeyBytes;
}
}
public static void main(String[] args) {
try {
String plaintext = "Hello, world!";
String secretKey = "ThisIsASecretKey123";
String encryptedText = encrypt(plaintext, secretKey);
System.out.println("Encrypted: " + encryptedText);
String decryptedText = decrypt(encryptedText, secretKey);
System.out.println("Decrypted: " + decryptedText);
} catch (Exception e) {
e.printStackTrace();
}
}
}






##### 5) Write Java/Python/C++ code to implement RSA Encryption Algorithm. 


#RSAExample.java




import java.security.KeyPair;
import java.security.KeyPairGenerator;
import java.security.PrivateKey;
import java.security.PublicKey;
import java.security.SecureRandom;
import java.security.Signature;
public class RSAExample {
public static KeyPair generateKeyPair() throws Exception {
SecureRandom secureRandom = new SecureRandom();
KeyPairGenerator keyPairGenerator = KeyPairGenerator.getInstance("RSA");
keyPairGenerator.initialize(2048, secureRandom);
return keyPairGenerator.generateKeyPair();
}
public static byte[] sign(PrivateKey privateKey, byte[] data) throws Exception {
Signature signature = Signature.getInstance("SHA256withRSA");
signature.initSign(privateKey);
signature.update(data);
return signature.sign();
}
public static boolean verify(PublicKey publicKey, byte[] data, byte[] signatureBytes) throws
Exception {
Signature signature = Signature.getInstance("SHA256withRSA");
signature.initVerify(publicKey);
signature.update(data);
return signature.verify(signatureBytes);
}
public static void main(String[] args) {
try {
// Generate key pair
KeyPair keyPair = generateKeyPair();
PublicKey publicKey = keyPair.getPublic();
PrivateKey privateKey = keyPair.getPrivate();
// Original data
String originalData = "Hello, world!";
byte[] originalDataBytes = originalData.getBytes();
// Sign the data with the private key
byte[] signatureBytes = sign(privateKey, originalDataBytes);
// Verify the signature using the public key
boolean isVerified = verify(publicKey, originalDataBytes, signatureBytes);
System.out.println("Original Data: " + originalData);
System.out.println("Signature Verified: " + isVerified);
} catch (Exception e) {
e.printStackTrace();
}
}
}



##### 6) Write code to Hellman Key Exchange mechanism and establish a secret key among 2 parties.




#Diffie-Hellman.html



<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta http-equiv="X-UA-Compatible" content="IE=edge">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Document</title>
</head>
<body>
<script>
// This program calculates the Key for two persons
// using the Diffie-Hellman Key exchange algorithm
// Power function to return value of a ^ b mod P
function power(a, b, p)
{
if (b == 1)
return a;
else
return((Math.pow(a, b)) % p);
}
// Driver code
var P, G, x, a, y, b, ka, kb;
// Both the persons will be agreed upon the
// public keys G and P
// A prime number P is taken
P = 23;
document.write("The value of P:" + P + "<br>");
// A primitive root for P, G is taken
G = 9;
document.write("The value of G:" + G + "<br>");
// Alice will choose the private key a
// a is the chosen private key
a = 4;
document.write("The private key a for Alice:" +
a + "<br>");
// Gets the generated key
x = power(G, a, P);
// Bob will choose the private key b
// b is the chosen private key
b = 3;
document.write("The private key b for Bob:" +
b + "<br>");
// Gets the generated key
y = power(G, b, P);
// Generating the secret key after the exchange
// of keys
ka = power(y, a, P); // Secret key for Alice
kb = power(x, b, P); // Secret key for Bob
document.write("Secret key for the Alice is:" +
ka + "<br>");
document.write("Secret key for the Bob is:" +
kb + "<br>");
</script>
</body>
</html>




#### 7) Write a program to Calculate the message digest of a text using the MD5 algorithm.  #####


# MD5Example.java



import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
public class MD5Example {
public static String calculateMD5(String text) throws NoSuchAlgorithmException {
// Create an instance of MessageDigest with MD5 algorithm
MessageDigest md = MessageDigest.getInstance("MD5");
// Convert the text to bytes and update the digest
md.update(text.getBytes());
// Get the digest bytes
byte[] digest = md.digest();
// Convert the digest bytes to a hexadecimal string representation
StringBuilder hexString = new StringBuilder();
for (byte b : digest) {
String hex = Integer.toHexString(0xFF & b);
if (hex.length() == 1) {
// Pad single digit hex values with leading zero
hexString.append('0');
}
hexString.append(hex);
}
return hexString.toString();
}
public static void main(String[] args) {
try {
String text = "Hello, world!";
String md5Digest = calculateMD5(text);
System.out.println("Text: " + text);
System.out.println("MD5 Digest: " + md5Digest);
} catch (NoSuchAlgorithmException e) {
e.printStackTrace();
}
}
}
