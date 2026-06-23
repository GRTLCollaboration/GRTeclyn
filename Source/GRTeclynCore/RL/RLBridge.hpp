#ifndef RL_BRIDGE_HPP_
#define RL_BRIDGE_HPP_

#include <AMReX_ParallelDescriptor.H>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>

#ifdef USE_RL
#include <zmq.hpp>
#endif

//! ZMQ binary bridge between GRTeclyn (REP) and Python RL agent (REQ).
class RLBridge
{
  public:
    static RLBridge &instance()
    {
        static RLBridge bridge;
        return bridge;
    }

    void configure(int port, int recv_timeout_ms, int send_timeout_ms)
    {
        m_port             = port;
        m_recv_timeout_ms  = recv_timeout_ms;
        m_send_timeout_ms  = send_timeout_ms;
    }

    void ensure_bound()
    {
#ifdef USE_RL
        if (m_bound)
        {
            return;
        }
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            m_ctx   = std::make_unique<zmq::context_t>(1);
            m_socket = std::make_unique<zmq::socket_t>(*m_ctx, zmq::socket_type::rep);
            m_socket->set(zmq::sockopt::rcvtimeo, m_recv_timeout_ms);
            m_socket->set(zmq::sockopt::sndtimeo, m_send_timeout_ms);
            const std::string addr = "tcp://127.0.0.1:" + std::to_string(m_port);
            m_socket->bind(addr);
            m_bound = true;
        }
#endif
    }

    std::vector<double> exchange(const std::vector<double> &obs, int action_dim,
                                 int timeout_ms)
    {
        (void)timeout_ms;
        std::vector<double> actions(static_cast<std::size_t>(action_dim), 0.0);

#ifdef USE_RL
        ensure_bound();
        if (amrex::ParallelDescriptor::IOProcessor() && m_socket)
        {
            if (!m_handshake_done)
            {
                zmq::message_t hello;
                const auto hello_result =
                    m_socket->recv(hello, zmq::recv_flags::none);
                if (!hello_result.has_value())
                {
                    return actions;
                }
                m_handshake_done = true;
            }

            zmq::message_t obs_msg(obs.size() * sizeof(double));
            std::memcpy(obs_msg.data(), obs.data(), obs.size() * sizeof(double));
            m_socket->send(obs_msg, zmq::send_flags::none);

            zmq::message_t reply;
            const auto recv_result =
                m_socket->recv(reply, zmq::recv_flags::none);
            if (recv_result.has_value() && reply.size() >=
                                                static_cast<std::size_t>(
                                                    action_dim) * sizeof(double))
            {
                std::memcpy(actions.data(), reply.data(),
                            static_cast<std::size_t>(action_dim) *
                                sizeof(double));
                if (std::isnan(actions[0]))
                {
                    m_terminate_requested = true;
                    actions.assign(static_cast<std::size_t>(action_dim), 0.0);
                }
            }
        }
#endif

#ifdef AMREX_USE_MPI
        amrex::ParallelDescriptor::Bcast(
            actions.data(), action_dim,
            amrex::ParallelDescriptor::IOProcessorNumber());
#endif
        return actions;
    }

    [[nodiscard]] bool terminate_requested() const
    {
        return m_terminate_requested;
    }

    void reset_terminate() { m_terminate_requested = false; }

  private:
    RLBridge() = default;

    int m_port{5555};
    int m_recv_timeout_ms{30000};
    int m_send_timeout_ms{30000};
    bool m_bound{false};
    bool m_handshake_done{false};
    bool m_terminate_requested{false};

#ifdef USE_RL
    std::unique_ptr<zmq::context_t> m_ctx;
    std::unique_ptr<zmq::socket_t> m_socket;
#endif
};

inline RLBridge &get_rl_bridge() { return RLBridge::instance(); }

#endif /* RL_BRIDGE_HPP_ */
