using UnityEngine;

public class CameraMovement : MonoBehaviour
{
    [SerializeField]
    private float moveSpeed = 10f;
    [SerializeField]
    private float lookSpeed = 2f;
    [SerializeField]
    private float verticalMoveSpeed = 5f;

    private float yaw = 0f;
    private float pitch = 0f;

    void Update()
    {
        float horizontal = Input.GetAxis("Horizontal"); 
        float vertical = Input.GetAxis("Vertical");     
        Vector3 direction = new Vector3(horizontal, 0, vertical);
        direction = transform.TransformDirection(direction);

        if (Input.GetKey(KeyCode.LeftControl))
        {
            transform.position += Vector3.down * verticalMoveSpeed * Time.deltaTime;
        }
        if (Input.GetKey(KeyCode.LeftShift))
        {
            transform.position += Vector3.up * verticalMoveSpeed * Time.deltaTime;
        }


        transform.position += direction * moveSpeed * Time.deltaTime;


        if (Input.GetMouseButton(1))
        {
            yaw += lookSpeed * Input.GetAxis("Mouse X");
            pitch -= lookSpeed * Input.GetAxis("Mouse Y");
            pitch = Mathf.Clamp(pitch, -90f, 90f);

            transform.eulerAngles = new Vector3(pitch, yaw, 0f);
        }
    }
}
